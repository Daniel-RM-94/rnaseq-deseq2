test_that("classify_de applies both padj and fold change thresholds", {
  status <- classify_de(
    log2fc = c(2, -2, 0.5, 3, NA, 2),
    padj = c(0.01, 0.01, 0.01, 0.2, 0.01, NA),
    alpha = 0.05, lfc_threshold = 1
  )
  expect_equal(as.character(status), c("Up", "Down", "NS", "NS", "NS", "NS"))
  expect_equal(levels(status), c("Up", "Down", "NS"))
})

test_that("keep_expressed_genes requires min_count in min_samples", {
  counts <- rbind(
    g1 = c(10, 10, 0, 0),
    g2 = c(10, 10, 10, 0),
    g3 = c(0, 0, 0, 500)
  )
  expect_equal(unname(keep_expressed_genes(counts, 10, 3)),
               c(FALSE, TRUE, FALSE))
})

test_that("prepare_sample_metadata sets the reference level", {
  samples <- data.frame(treatment = c("trt", "untrt", "trt", "untrt"),
                        cell = c("a", "a", "b", "b"))
  design <- list(formula = "~ cell + treatment", factor = "treatment",
                 reference = "untrt", treatment = "trt")

  out <- prepare_sample_metadata(samples, design)
  expect_equal(levels(out$treatment), c("untrt", "trt"))
  expect_s3_class(out$cell, "factor")
})

test_that("prepare_sample_metadata validates the design", {
  samples <- data.frame(treatment = c("trt", "untrt"))
  design <- list(formula = "~ cell + treatment", factor = "treatment",
                 reference = "untrt", treatment = "trt")
  expect_error(prepare_sample_metadata(samples, design), "not found")

  design$formula <- "~ treatment"
  design$reference <- "control"
  expect_error(prepare_sample_metadata(samples, design), "Levels not found")
})

test_that("strip_ensembl_version removes version suffixes only", {
  expect_equal(strip_ensembl_version(c("ENSG00000000003.15", "ENSG00000000005")),
               c("ENSG00000000003", "ENSG00000000005"))
})

test_that("DESeq2 results are formatted into a tidy, sorted table", {
  skip_if_not_installed("DESeq2")
  set.seed(1)
  dds <- DESeq2::makeExampleDESeqDataSet(n = 500, m = 6)
  dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
  res <- DESeq2::results(dds, contrast = c("condition", "B", "A"))

  de_table <- format_de_results(res, res, annotation = NULL,
                                alpha = 0.1, lfc_threshold = 0)
  expect_s3_class(de_table, "tbl_df")
  expect_equal(nrow(de_table), 500)
  expect_true(all(c("gene_id", "log2FoldChange", "stat", "padj", "status")
                  %in% names(de_table)))
  expect_false(is.unsorted(de_table$padj, na.rm = TRUE))

  de_summary <- summarise_de(de_table)
  expect_equal(de_summary$total_significant, sum(de_table$status != "NS"))

  gene_list <- ranked_gene_list(de_table, keytype = "SYMBOL")
  expect_false(is.unsorted(rev(gene_list)))
  expect_false(anyDuplicated(names(gene_list)) > 0)
})
