# Register global variables to avoid "no visible binding for global variable"
# NOTES from lintr / R CMD check. Add any unquoted column names used in
# NSE (dplyr, data.table, ggplot2 aes, etc.)
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      # gene / DE result columns used bare in dplyr / ggplot2 NSE
      "Geneid", "symbol", "padj", "pval", "log2FC", "log2FC_shrunk",
      "regulated", "DE", "FC", "stat", "baseMean",
      
      # contrast / pathway columns used bare in NSE
      "contrast", "pathway",
      
      # GSEA result columns used bare in NSE
      "NES", "ES", "size", "leadingEdge", "pctVar", "log2err",
      
      # genomic / neighbourhood columns used bare in NSE
      "start", "end", "strand", "track",
      
      # ggplot2 tidy-eval helper
      ".data"
    )
  )
}

