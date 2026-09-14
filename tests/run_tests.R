#!/usr/bin/env Rscript
# Ejecuta los tests unitarios y termina con error si alguno falla.

testthat::test_dir(
  here::here("tests", "testthat"),
  reporter = testthat::ProgressReporter$new(show_praise = FALSE),
  stop_on_failure = TRUE
)
