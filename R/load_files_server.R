load_files_server <- function(input, output, session, state) {

  # ---- Helpers ----------------------------------------------------------------

  # Locate a built-in annotation file for the given organism name.
  find_annotation_path <- function(organism) {
    # Try installed-package path first.
    path <- system.file("annotations",
                        paste0(organism, "_annot.rds"),
                        package = "RNASeqApp")
    if (nchar(path) > 0) return(path)

    # Fall back to relative path when running the app directly from source.
    local_path <- file.path("inst", "annotations",
                            paste0(organism, "_annot.rds"))
    if (file.exists(local_path)) return(local_path)

    return(NULL)
  }

  # Extract all contrasts from SE rowData into a long-format data frame.
  extract_contrasts <- function(se) {
    rd       <- as.data.frame(SummarizedExperiment::rowData(se))
    gene_ids <- rownames(se)

    lfc_cols <- grep("^log2FoldChange_", colnames(rd), value = TRUE)
    if (length(lfc_cols) == 0) {
      stop("No contrast columns found in rowData. ",
           "Expected columns like 'log2FoldChange_ContrastName'.")
    }

    contrast_names <- sub("^log2FoldChange_", "", lfc_cols)

    de_list <- lapply(contrast_names, function(cname) {
      get_col <- function(prefix) {
        col <- paste0(prefix, cname)
        if (col %in% colnames(rd)) rd[[col]] else rep(NA_real_, nrow(rd))
      }
      data.frame(
        Geneid        = gene_ids,
        contrast      = cname,
        baseMean      = get_col("baseMean_"),
        log2FC        = get_col("log2FoldChange_"),
        log2FC_shrunk = get_col("log2FoldChange_shrunk_"),
        lfcSE         = get_col("lfcSE_"),
        stat          = get_col("stat_"),
        pvalue        = get_col("pvalue_"),
        padj          = get_col("padj_"),
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, de_list)
  }

  # Merge an annotation data frame into the contrast data frame.
  merge_annotation <- function(de_df, annot) {
    if (!"Geneid" %in% colnames(annot)) {
      showModal(modalDialog(
        title = "Annotation error",
        "Annotation file is missing required 'Geneid' column.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(de_df)
    }

    # Store full annotation for export functionality.
    state$annotation_df(annot)

    dplyr::left_join(de_df, annot, by = "Geneid")
  }

  # ---- SE file upload ---------------------------------------------------------

  observeEvent(input$seFile, {
    req(input$seFile)

    se <- tryCatch(
      readRDS(input$seFile$datapath),
      error = function(e) {
        showModal(modalDialog(
          title = "File error",
          paste("Error reading file:", e$message),
          easyClose = TRUE,
          footer = NULL
        ))
        NULL
      }
    )
    req(se)

    if (!inherits(se, "SummarizedExperiment")) {
      showModal(modalDialog(
        title = "File error",
        "Uploaded file is not a SummarizedExperiment object.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }

    # Store SE (compatible with DESeq2/SummarizedExperiment accessor functions).
    state$dds_obj(se)

    # Extract contrasts from rowData.
    de_df <- tryCatch(
      extract_contrasts(se),
      error = function(e) {
        showModal(modalDialog(
          title = "Contrast extraction error",
          e$message,
          easyClose = TRUE,
          footer = NULL
        ))
        NULL
      }
    )
    req(de_df)

    # Extract variance partition if varpart_* columns are present.
    rd <- as.data.frame(SummarizedExperiment::rowData(se))
    varpart_cols <- grep("^varpart_", colnames(rd), value = TRUE)
    if (length(varpart_cols) > 0) {
      vp_mat <- as.data.frame(rd[, varpart_cols, drop = FALSE])
      colnames(vp_mat) <- sub("^varpart_", "", colnames(vp_mat))
      rownames(vp_mat) <- rownames(se)
      state$varpart_obj(list(varPart = vp_mat))
    }

    # Determine organism: prefer SE metadata, fall back to text input.
    se_organism <- tryCatch(
      metadata(se)$organism,
      error = function(e) NULL
    )
    organism <- if (!is.null(se_organism) && nzchar(trimws(se_organism))) {
      se_organism
    } else {
      input$organism_name
    }

    if (is.null(organism) || !nzchar(trimws(organism))) {
      showNotification(
        "No organism specified. Enter an organism name to load annotations.",
        type = "warning", duration = 8
      )
      state$de_df(de_df)
      output$annotationStatus <- renderUI({
        tagList(
          tags$p(icon("exclamation-triangle"), "No organism specified.",
                 style = "color: orange; padding-left: 15px;"),
          fileInput("annotFile", "Upload Annotation (.rds)", accept = ".rds")
        )
      })
      return(NULL)
    }

    # Try to load built-in annotation.
    annot_path <- find_annotation_path(organism)

    if (!is.null(annot_path)) {
      annot <- tryCatch(
        readRDS(annot_path),
        error = function(e) {
          showNotification(paste("Annotation load error:", e$message),
                           type = "error")
          NULL
        }
      )
      if (!is.null(annot)) {
        de_df <- merge_annotation(de_df, annot)
      }
      state$de_df(de_df)
      output$annotationStatus <- renderUI({
        tags$p(icon("check-circle"),
               paste("Annotation loaded:", basename(annot_path)),
               style = "color: green; padding-left: 15px;")
      })
    } else {
      # No built-in annotation found -- show upload UI.
      state$de_df(de_df)
      output$annotationStatus <- renderUI({
        tagList(
          tags$p(icon("exclamation-triangle"),
                 paste("No built-in annotation found for:", organism),
                 style = "color: orange; padding-left: 15px;"),
          fileInput("annotFile", "Upload Annotation (.rds)", accept = ".rds")
        )
      })
    }
  })

  # ---- Custom annotation upload -----------------------------------------------

  observeEvent(input$annotFile, {
    req(input$annotFile, state$de_df())

    annot <- tryCatch(
      readRDS(input$annotFile$datapath),
      error = function(e) {
        showModal(modalDialog(
          title = "Annotation error",
          paste("Error reading annotation file:", e$message),
          easyClose = TRUE,
          footer = NULL
        ))
        NULL
      }
    )
    req(annot)

    de_df <- merge_annotation(state$de_df(), annot)
    state$de_df(de_df)
    output$annotationStatus <- renderUI({
      tags$p(icon("check-circle"), "Custom annotation loaded successfully.",
             style = "color: green; padding-left: 15px;")
    })
    showNotification("Custom annotation loaded successfully.", type = "message")
  })
}