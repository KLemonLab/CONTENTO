load_files_server <- function(input, output, session, state) {
  
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
    
    
    required_cols <- c("contrast", "Geneid", "symbol", "padj",
                       "log2FoldChange", "log2FoldChange_shrunk")
    
    if (!all(required_cols %in% colnames(df))) {
      showModal(modalDialog(
        title = "File error",
        paste0("File missing required columns: ", paste(setdiff(required_cols, colnames(df)), collapse = ", ")),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }
    
    # Rename columns
    if ("log2FoldChange" %in% colnames(df)) {
      df <- dplyr::rename(df, log2FC = "log2FoldChange")
    }
    if ("log2FoldChange_shrunk" %in% colnames(df)) {
      df <- dplyr::rename(df, log2FC_shrunk = "log2FoldChange_shrunk")
    }
    
    state$de_df(df)
  })
  
  # VarPart (.rds)
  observeEvent(input$varPartFile, {
    req(input$varPartFile)
    state$varpart_obj(readRDS(input$varPartFile$datapath))
  })
}