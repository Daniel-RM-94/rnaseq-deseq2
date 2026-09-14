#!/usr/bin/env Rscript
# Ejecuta los pasos del análisis en orden, cada uno en una sesión de R nueva,
# genera el informe HTML y guarda las versiones del software usado.
#
# Uso:
#   Rscript scripts/00_run_pipeline.R              # análisis e informe
#   Rscript scripts/00_run_pipeline.R --no-report  # solo el análisis

args <- commandArgs(trailingOnly = TRUE)
render_report <- !"--no-report" %in% args

source(here::here("R", "utils.R"))
cfg <- load_config()

steps <- c(
  "01_prepare_data.R",
  "02_quality_control.R",
  "03_differential_expression.R",
  "04_functional_enrichment.R"
)

rscript <- file.path(R.home("bin"), "Rscript")
pipeline_start <- Sys.time()

for (step in steps) {
  log_step("==> Running ", step)
  step_start <- Sys.time()
  status <- system2(rscript, shQuote(here::here("scripts", step)))
  if (!identical(as.integer(status), 0L)) {
    stop("Pipeline step failed: ", step, " (exit status ", status, ")",
         call. = FALSE)
  }
  log_step(sprintf("<== %s finished in %.1f s", step,
                   as.numeric(difftime(Sys.time(), step_start, units = "secs"))))
}

if (render_report) {
  log_step("==> Rendering HTML report")
  rmarkdown::render(
    here::here("reports", "rnaseq_report.Rmd"),
    output_dir = ensure_dir(cfg$paths$reports),
    knit_root_dir = here::here(),
    intermediates_dir = tempdir(),
    envir = new.env(),
    quiet = TRUE
  )
}

# Copia las figuras principales a docs/ para mostrarlas en el README.
showcase <- c("qc_pca.png", "de_volcano.png", "de_heatmap_top_genes.png",
              "go_gsea_dotplot.png")
showcase_src <- file.path(here::here(cfg$paths$figures), showcase)
showcase_src <- showcase_src[file.exists(showcase_src)]
invisible(file.copy(showcase_src, ensure_dir(cfg$paths$showcase_figures),
                    overwrite = TRUE))

session_file <- file.path(ensure_dir(cfg$paths$results), "session_info.txt")
writeLines(
  utils::capture.output(sessioninfo::session_info(
    pkgs = c("DESeq2", "apeglm", "clusterProfiler", "enrichplot",
             "org.Hs.eg.db", "airway", "ggplot2", "pheatmap"),
    dependencies = FALSE
  )),
  session_file
)

log_step(sprintf("Pipeline completed in %.1f min",
                 as.numeric(difftime(Sys.time(), pipeline_start,
                                     units = "mins"))))
