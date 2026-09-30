# TEST ANNOTATIONS

This directory contains intentionally malformed and edge-case annotation files to use for testing CONTENTO's annotation validation and GSEA behavior.

Files:
- annotation_good.csv: Well-formed annotation with `Geneid` (character) and `Regulator` using the expected `!!!` separator for multi-value entries.
- annotation_missing_geneid.csv: Missing the required `Geneid` column (tests error message shown in CONTENTO when joining annotations).
- annotation_geneid_numeric.csv: `Geneid` column present but stored as numeric values (tests type-check failure in merge_annotation).
- annotation_multiline_names.csv: `Regulator` column contains multiple gene names separated by literal newlines instead of the `!!!` separator (tests parsing/cleanup routines).
- annotation_zerosets.csv: Annotation column `tiny_group` contains unique values for each gene so selecting it for annotation-based GSEA should yield zero gene sets after standard size filtering (min size >= 2).

R script:
- generate_annotations.R: Convenience script to convert the CSVs into `.rds` files under `inst/extdata/test_annotations/rds/` and to force specific column types (e.g., numeric Geneid).

Usage:
1. From the package root, run:
   Rscript inst/extdata/test_annotations/generate_annotations.R

2. Inspect resulting `.rds` files in `inst/extdata/test_annotations/rds/` and load them into CONTENTO (or use within unit tests) to verify validation and behavior.
