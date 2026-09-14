#!/usr/bin/env Rscript
# Ajusta el modelo binomial negativo de DESeq2, evalúa el contraste configurado
# (test de Wald), contrae los log2 fold change y anota los resultados.
#
# Entrada: data/processed/dds_filtered.rds
# Salida:  data/processed/{dds_deseq,vsd_design,de_results}.rds
#          results/tables/de_results_{all,significant}.tsv, de_summary.tsv
#          results/figures/de_*.png

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
design_cfg <- cfg$design
de_cfg <- cfg$de
plots_cfg <- cfg$plots

dds <- readRDS(file.path(processed_dir, "dds_filtered.rds"))

# 1. Ajuste del modelo ----
log_step("Running DESeq(): size factors, dispersions and GLM fit")
dds <- DESeq(dds, quiet = TRUE)
saveRDS(dds, file.path(processed_dir, "dds_deseq.rds"))

# 2. Test del contraste de interés ----
contrast <- c(design_cfg$factor, design_cfg$treatment, design_cfg$reference)
log_step("Testing contrast: ", design_cfg$treatment, " vs ",
         design_cfg$reference, " (", design_cfg$factor, ")")
res <- results(dds, contrast = contrast, alpha = de_cfg$alpha)

# 3. Contracción del log2 fold change ----
res_shrunk <- switch(
  de_cfg$shrinkage,
  apeglm = {
    coef_name <- paste(design_cfg$factor, design_cfg$treatment, "vs",
                       design_cfg$reference, sep = "_")
    if (!coef_name %in% resultsNames(dds)) {
      stop("Coefficient '", coef_name, "' not found. Available: ",
           paste(resultsNames(dds), collapse = ", "), call. = FALSE)
    }
    lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
  },
  normal = lfcShrink(dds, contrast = contrast, type = "normal", quiet = TRUE),
  none = res
)

# 4. Anotación y tablas de resultados ----
orgdb <- load_orgdb(cfg$annotation$orgdb)
annotation <- annotate_genes(rownames(res), orgdb, cfg$annotation$keytype)

de_table <- format_de_results(res, res_shrunk, annotation,
                              alpha = de_cfg$alpha,
                              lfc_threshold = de_cfg$lfc_threshold)
de_summary <- summarise_de(de_table)
log_step(sprintf("Significant genes: %d up, %d down (padj < %s, |log2FC| >= %s)",
                 de_summary$up, de_summary$down, de_cfg$alpha,
                 de_cfg$lfc_threshold))

saveRDS(de_table, file.path(processed_dir, "de_results.rds"))
write_table(de_table, tables_dir, "de_results_all.tsv")
write_table(dplyr::filter(de_table, status != "NS"), tables_dir,
            "de_results_significant.tsv")
write_table(de_summary, tables_dir, "de_summary.tsv")

# 5. Figuras ----
# VST que tiene en cuenta el diseño, para visualizar los genes DE.
vsd <- variance_stabilize(dds, blind = FALSE)
saveRDS(vsd, file.path(processed_dir, "vsd_design.rds"))

save_base_plot(function() plotDispEsts(dds, main = "Estimación de dispersión"),
               figures_dir, "de_dispersion.png", dpi = plots_cfg$dpi)

save_plot(plot_ma(de_table, de_cfg$lfc_threshold),
          figures_dir, "de_ma_plot.png", dpi = plots_cfg$dpi)

save_plot(plot_volcano(de_table, de_cfg$alpha, de_cfg$lfc_threshold,
                       n_labels = plots_cfg$n_volcano_labels),
          figures_dir, "de_volcano.png", width = 8, height = 6,
          dpi = plots_cfg$dpi)

heatmap_file <- file.path(figures_dir, "de_heatmap_top_genes.png")
if (!is.null(plot_top_genes_heatmap(vsd, de_table, plots_cfg$n_top_heatmap,
                                    plots_cfg$annotation_cols,
                                    filename = heatmap_file))) {
  log_step("Saved figure: ", heatmap_file)
}

save_plot(plot_gene_counts(dds, de_table, design_cfg$factor,
                           plots_cfg$shape_by, n = plots_cfg$n_count_plots),
          figures_dir, "de_top_genes_counts.png", width = 9, height = 6,
          dpi = plots_cfg$dpi)
