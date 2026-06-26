# ==============================================================================
# DEV PRELOAD — mirrors load_files_server SE upload logic without fileInput.
# Only active when DEV_CONFIG$enabled is TRUE.
# ==============================================================================

dev_preload_server <- function(input, output, session, state, dev_config) {

  observe({
    req(dev_config$enabled, nzchar(dev_config$se_path %||% ""))

    se_path <- dev_config$se_path

    if (!file.exists(se_path)) {
      message("[DEV] SE file not found: ", se_path)
      return()
    }

    message("[DEV] Loading SE from: ", se_path)
    se <- tryCatch(readRDS(se_path), error = function(e) {
      message("[DEV] Failed to read SE: ", e$message); NULL
    })
    if (is.null(se) || !inherits(se, "SummarizedExperiment")) {
      message("[DEV] Object is not a SummarizedExperiment — aborting preload.")
      return()
    }

    state$se_obj(se)

    # ── DE results ────────────────────────────────────────────────────────────
    de_df <- tryCatch(extract_de_results(se), error = function(e) {
      message("[DEV] extract_de_results failed: ", e$message); NULL
    })
    if (is.null(de_df)) return()

    # ── Variance partition ────────────────────────────────────────────────────
    rd           <- as.data.frame(rowData(se))
    varpart_cols <- grep("^varpart_", colnames(rd), value = TRUE)
    if (length(varpart_cols) > 0) {
      vp_mat           <- as.data.frame(rd[, varpart_cols, drop = FALSE])
      colnames(vp_mat) <- sub("^varpart_", "", colnames(vp_mat))
      rownames(vp_mat) <- rownames(se)
      state$varpart_obj(list(varPart = vp_mat))
    }

    # ── Organism ──────────────────────────────────────────────────────────────
    se_organism <- tryCatch(metadata(se)$organism, error = function(e) NULL)
    organism    <- if (!is.null(se_organism) && nzchar(trimws(se_organism))) se_organism else NULL
    state$se_organism(organism)

    # ── Annotation ────────────────────────────────────────────────────────────
    se_annotation <- tryCatch(metadata(se)$annotation, error = function(e) NULL)
    annotation    <- if (!is.null(se_annotation) && nzchar(trimws(se_annotation))) se_annotation else NULL

    annot <- NULL

    if (!is.null(annotation)) {
      annot_path <- if (!is.null(dev_config$annotation_path) && nzchar(dev_config$annotation_path)) {
        dev_config$annotation_path
      } else {
        find_annotation_path(annotation)
      }

      if (!is.null(annot_path) && file.exists(annot_path)) {
        annot <- tryCatch(readRDS(annot_path), error = function(e) {
          message("[DEV] Failed to read annotation: ", e$message); NULL
        })

        if (!is.null(annot)) {
          merged_result <- merge_annotation(de_df, annot)
          if (!is.null(merged_result$result)) {
            merged      <- merged_result$result
            annot_cols  <- colnames(annot)
            sym_primary <- default_symbol_col(annot_cols)
            de_df <- apply_symbol(merged, sym_primary, "Geneid")
            annot <- apply_symbol(annot,  sym_primary, "Geneid")
            state$annotation_df(annot)
            message("[DEV] Annotation merged. Symbol column: ", sym_primary)
          } else {
            message("[DEV] Annotation merge failed: ", merged_result$message)
          }
        }
      } else {
        message("[DEV] Annotation path not found for: ", annotation)
      }
    }

    state$de_df(de_df)

    # ── GSEA / MSigDB derived state ───────────────────────────────────────────
    if (!is.null(state$annotation_df())) {
      state$available_gsea_columns(get_gsea_columns(state$annotation_df()))
      state$ensembl_col(detect_ensembl_col(state$annotation_df()))
    }
    db_override <- tryCatch(metadata(se)$db_species, error = function(e) NULL)
    state$db_species(get_msigdbr_db_species(organism, db_override))

    message("[DEV] State populated. Contrasts available: ",
            paste(unique(de_df$contrast), collapse = ", "))

    # ── Pre-select default contrasts ─────────────────────────────────────────
    # contrastSelect and multiContrastCheckboxes are renderUI blocks that only
    # render when their tab is *visited* — so input$contrast and
    # input$compare_contrasts do not exist at startup.
    # We store the desired selections in a reactiveVal and apply them lazily
    # whenever each input first appears (tab visit).
    if (length(dev_config$default_contrasts) > 0) {
      valid <- intersect(dev_config$default_contrasts, unique(de_df$contrast))
      if (length(valid) == 0) {
        message("[DEV] None of default_contrasts matched: ",
                paste(dev_config$default_contrasts, collapse = ", "),
                "\n  Available: ", paste(unique(de_df$contrast), collapse = ", "))
      } else {
        all_contrasts <- unique(de_df$contrast)
        message("[DEV] Will pre-select: ", paste(valid, collapse = ", "))

        # ── single-contrast selectInput (Explore by Contrast tab) ─────────────
        single_done <- FALSE
        observe({
          req(input$contrast)
          if (single_done) return()
          single_done <<- TRUE
          updateSelectInput(session, "contrast", selected = valid[1])
          message("[DEV] contrast selectInput updated to: ", valid[1])
        })

        # ── multi-contrast checkboxes (Compare tab) ────────────────────────────
        # multiContrastCheckboxes renderUI uses isolate(input$compare_contrasts)
        # to preserve selection on re-render. We need to fire AFTER that renderUI
        # has run (i.e. after the Compare tab is first visited).
        # Watching input$tabs for "compare" is the reliable trigger.
        multi_done <- FALSE
        observe({
          req(input$tabs == "compare")
          if (multi_done) return()
          # Give the renderUI one extra flush to finish drawing the checkboxes
          invalidateLater(200, session)
          isolate({
            if (multi_done) return()
            multi_done <<- TRUE
            updateCheckboxGroupInput(session, "compare_contrasts",
                                     choices  = all_contrasts,
                                     selected = valid)
            message("[DEV] compare_contrasts updated to: ", paste(valid, collapse = ", "))
          })
        })
      }
    }

    # ── Optionally jump straight to a tab ─────────────────────────────────────
    if (!is.null(dev_config$start_tab) && nzchar(dev_config$start_tab)) {
      session$onFlushed(function() {
        updateTabItems(session, "tabs", selected = dev_config$start_tab)
      }, once = TRUE)
    }

  })
}
