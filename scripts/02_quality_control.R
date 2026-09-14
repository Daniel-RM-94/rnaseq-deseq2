#!/usr/bin/env Rscript
# Control de calidad: tamaño de librería, distribución de conteos, PCA y
# distancias entre muestras.
#
# Entrada: data/processed/dds_filtered.rds
# Salida:  data/processed/vsd_blind.rds
#          results/tables/qc_sample_metrics.tsv
#          results/figures/qc_*.png

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
figures_dir <- ensure_dir(cfg$paths$figures)
plots_cfg <- cfg$plots

dds <- readRDS(file.path(processed_dir, "dds_filtered.rds"))
dds <- estimateSizeFactors(dds)

# 1. Métricas por muestra ----
qc_metrics <- tibble::tibble(
  sample_id = colnames(dds),
  total_counts = colSums(counts(dds)),
  detected_genes = colSums(counts(dds) > 0),
  size_factor = sizeFactors(dds)
)
write_table(qc_metrics, tables_dir, "qc_sample_metrics.tsv")

# 2. Transformación VST sin usar el diseño (blind), como corresponde al QC ----
log_step("Applying variance-stabilising transformation (blind = TRUE)")
vsd <- variance_stabilize(dds, blind = TRUE)
saveRDS(vsd, file.path(processed_dir, "vsd_blind.rds"))

# 3. Figuras ----
save_plot(plot_library_sizes(dds, plots_cfg$color_by),
          figures_dir, "qc_library_sizes.png", dpi = plots_cfg$dpi)

save_plot(plot_count_distribution(dds, plots_cfg$color_by),
          figures_dir, "qc_count_distribution.png",
          width = 8, height = 7, dpi = plots_cfg$dpi)

save_plot(plot_pca(vsd, plots_cfg$color_by, plots_cfg$shape_by),
          figures_dir, "qc_pca.png", width = 7, height = 5,
          dpi = plots_cfg$dpi)

plot_sample_distances(vsd, plots_cfg$annotation_cols,
                      filename = file.path(figures_dir,
                                           "qc_sample_distances.png"))
log_step("Saved figure: ", file.path(figures_dir, "qc_sample_distances.png"))
