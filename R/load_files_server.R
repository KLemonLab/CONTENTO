load_files_server <- function(input, output, session, state) {
  
  # Required columns
  required_cols_bacteria <- c("contrast", "Geneid", "symbol", "padj",
                              "log2FoldChange", "log2FoldChange_shrunk", 
                              "start", "end", "strand", "biotype")
  required_cols_human <- c("contrast", "Geneid", "symbol", "padj",
                           "log2FoldChange", "log2FoldChange_shrunk")
  
  # DDS (.rds)
  observeEvent(input$ddsFile, {
    req(input$ddsFile)
    obj <- readRDS(input$ddsFile$datapath)
    if (!("DESeqDataSet" %in% class(obj))) {
      showModal(modalDialog(
        title = "File error",
        "Uploaded DDS file is not a DESeqDataSet object.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }
    state$dds_obj(obj)
  })
  
  # DE contrasts (.rds)
  observeEvent(input$deFile, {
    req(input$deFile)
    df <- readRDS(input$deFile$datapath)
    
    # Select required columns based on organism
    cols_to_check <- if (input$organism == "Bacteria") {
      required_cols_bacteria
    } else {
      required_cols_human
    }
    
    # Check columns
    if (!all(cols_to_check %in% colnames(df))) {
      showModal(modalDialog(
        title = "File error",
        paste0("File missing required columns: ", 
               paste(setdiff(cols_to_check, colnames(df)), collapse = ", ")),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }
    
    # Rename columns if needed
    if ("log2FoldChange" %in% colnames(df)) {
      df <- dplyr::rename(df, log2FC = "log2FoldChange")
    }
    if ("log2FoldChange_shrunk" %in% colnames(df)) {
      df <- dplyr::rename(df, log2FC_shrunk = "log2FoldChange_shrunk")
    }
    
    state$de_df(df)
    
    # Create annotation dataframe
    deseq_cols <- c("contrast", "baseMean", "log2FC", "lfcSE", "stat", 
                    "pvalue", "padj", "log2FC_shrunk", "sign", "DE", "regulated")
    
    annotation_cols <- setdiff(names(df), deseq_cols)
    
    annotation_df <- df %>%
      select(all_of(annotation_cols)) %>%
      distinct()
    
    state$annotation_df(annotation_df)
  })
  
  # VarPart (.rds)
  observeEvent(input$varPartFile, {
    req(input$varPartFile)
    state$varpart_obj(readRDS(input$varPartFile$datapath))
  })
}