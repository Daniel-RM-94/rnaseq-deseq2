#!/usr/bin/env Rscript
# Importa los conteos y los metadatos, crea el DESeqDataSet y filtra los genes
# poco expresados.
#
# Entrada: dataset 'airway' o archivos TSV propios (sección input de la config)
# Salida:  data/processed/dds_filtered.rds
#          results/tables/sample_metadata.tsv, filtering_summary.tsv

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
})

source(here::here("R", "utils.R"))
load_project_functions()

cfg <- load_config()
set.seed(cfg$project$seed)
processed_dir <- ensure_dir(cfg$paths$processed)
tables_dir <- ensure_dir(cfg$paths$tables)

# 1. Importar conteos y metadatos ----
if (cfg$input$source == "airway") {
  log_step("Loading the 'airway' dataset (Himes et al. 2014, GSE52778)")
  data("airway", package = "airway", envir = environment())
  counts_mat <- assay(airway, "counts")
  samples <- as.data.frame(colData(airway))
  samples <- cbind(sample_id = rownames(samples), samples)

  if (isTRUE(cfg$input$export_example_files)) {
    raw_dir <- ensure_dir("data/raw")
    readr::write_tsv(tibble::as_tibble(counts_mat, rownames = "gene_id"),
                     file.path(raw_dir, "airway_counts.tsv.gz"))
    readr::write_tsv(samples, file.path(raw_dir, "airway_samples.tsv"))
    log_step("Example input files exported to ", raw_dir)
  }
} else {
  log_step("Reading user-provided count matrix and sample sheet")
  counts_mat <- read_count_matrix(here::here(cfg$input$counts_file))
  samples <- read_sample_sheet(here::here(cfg$input$samples_file),
                               cfg$input$sample_id_column)
}

# 2. Validar y crear el DESeqDataSet ----
validate_count_inputs(counts_mat, samples)
samples <- samples[colnames(counts_mat), , drop = FALSE]
samples <- prepare_sample_metadata(samples, cfg$design)
storage.mode(counts_mat) <- "integer"

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = samples,
  design = stats::as.formula(cfg$design$formula)
)
n_imported <- nrow(dds)
log_step(sprintf("Imported %d genes x %d samples", n_imported, ncol(dds)))

# 3. Filtrar genes poco expresados ----
keep <- keep_expressed_genes(counts(dds),
                             min_count = cfg$filtering$min_count,
                             min_samples = cfg$filtering$min_samples)
dds <- dds[keep, ]
log_step(sprintf("Retained %d of %d genes (>= %d counts in >= %d samples)",
                 nrow(dds), n_imported, cfg$filtering$min_count,
                 cfg$filtering$min_samples))

# 4. Guardar resultados ----
saveRDS(dds, file.path(processed_dir, "dds_filtered.rds"))
write_table(samples, tables_dir, "sample_metadata.tsv")
write_table(
  tibble::tibble(
    stage = c("imported", "retained_after_filtering"),
    n_genes = c(n_imported, nrow(dds))
  ),
  tables_dir, "filtering_summary.tsv"
)
