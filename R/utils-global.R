# Register global variables to avoid "no visible binding for global variable"
# NOTES from lintr / R CMD check. Add any unquoted column names used in
# NSE (dplyr, data.table, ggplot2 aes, etc.)
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      # gene / DE result columns
      "Geneid", "symbol", "gene", "padj", "pval", "log2FC", "log2FC_shrunk",
      "regulated", "DE",
      
      # contrast / pathway / GSEA names
      "contrast", "contrast_str", "file_base", "pathway_str", "pathways_list",
      "pathway",  "vst_matrix", "gesecaRes", "tableplot",
      
      # matrices / leading edge / upset
      "leadingEdge", "leading_edge_matrix", "leading_edge_data", "upset_df",
      "padj_mat", "nes_mat", "n_pathways",
      
      # formatting / display related
      "pctVar", "log2err",
      
      # plotting / selection helpers
      "x_col", "color_col", "shape_col",
      
      ".data"
      
    )
  )
}

# source("R/utils-global.R")
# lintr:::addin_lint()
