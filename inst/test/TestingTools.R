### TEST ###

library(here)
source(here("R/utils.R"))

Sau_KPL4403_annot <- readRDS("inst/annotations/Sau_KPL4403_annot.rds")
SE_5772AIB_LAC <- readRDS("inst/test/5772AIB_LAC_SE.rds")
contrasts = metadata(SE_5772AIB_LAC)$contrasts
DE = extract_de_results(SE_5772AIB_LAC)
DE <- apply_symbol(DE, "Geneid", "Geneid")
selected_data <- add_de_flags(DE, lfc_cut = 1, padj_cut = 0.05)
compare_data <- build_compare_table(selected_data, contrasts)

Human_annot <- readRDS("inst/annotations/Homo_sapiens.GRCh38.113.rds")
SE_5772AIB_LAC <- readRDS("inst/test/RSVBac_human_SE.rds")
