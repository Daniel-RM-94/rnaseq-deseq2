test_that("project configuration loads and is consistent", {
  cfg <- load_config()
  expect_true(cfg$input$source %in% c("airway", "files"))
  expect_true(cfg$design$factor %in% design_variables(cfg$design$formula))
  expect_false(identical(cfg$design$reference, cfg$design$treatment))
  expect_gt(cfg$de$alpha, 0)
  expect_lt(cfg$de$alpha, 1)
})

test_that("load_config fails clearly when the file does not exist", {
  expect_error(load_config(tempfile(fileext = ".yml")), "not found")
})

test_that("validate_count_inputs accepts well-formed inputs", {
  counts <- matrix(c(0, 5, 10, 20), nrow = 2,
                   dimnames = list(c("g1", "g2"), c("s1", "s2")))
  samples <- data.frame(group = c("a", "b"), row.names = c("s2", "s1"))
  expect_true(validate_count_inputs(counts, samples))
})

test_that("validate_count_inputs rejects invalid counts", {
  samples <- data.frame(group = c("a", "b"), row.names = c("s1", "s2"))
  make_counts <- function(values, samples = c("s1", "s2")) {
    matrix(values, nrow = 2, dimnames = list(c("g1", "g2"), samples))
  }

  expect_error(validate_count_inputs(make_counts(c(-1, 2, 3, 4)), samples),
               "negative")
  expect_error(validate_count_inputs(make_counts(c(1.5, 2, 3, 4)), samples),
               "integer")
  expect_error(validate_count_inputs(make_counts(c(NA, 2, 3, 4)), samples),
               "missing")
  expect_error(validate_count_inputs(make_counts(1:4, c("s1", "s3")), samples),
               "do not match")
})

test_that("read_count_matrix uses the first column as gene identifiers", {
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path))
  writeLines(c("gene_id\ts1\ts2", "g1\t1\t2", "g2\t3\t4"), path)

  counts <- read_count_matrix(path)
  expect_equal(dim(counts), c(2L, 2L))
  expect_equal(rownames(counts), c("g1", "g2"))
  expect_equal(colnames(counts), c("s1", "s2"))
  expect_equal(unname(counts["g2", "s2"]), 4)
})

test_that("read_count_matrix rejects duplicated gene identifiers", {
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path))
  writeLines(c("gene_id\ts1", "g1\t1", "g1\t3"), path)
  expect_error(read_count_matrix(path), "Duplicated")
})
