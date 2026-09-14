#!/usr/bin/env Rscript
# Interpretación funcional con Gene Ontology:
#   - ORA: términos GO sobrerrepresentados en los genes sobre y subexpresados
#   - GSEA: enriquecimiento en el ranking completo de genes (estadístico de Wald)
#
# Entrada: data/processed/de_results.rds
# Salida:  data/processed/enrichment.rds
#          results/tables/go_{ora_up,ora_down,gsea}.tsv
#          results/figures/go_*.png

source(here::here("R", "utils.R"))
cfg <- load_config()

if (!isTRUE(cfg$enrichment$run)) {
  log_step("Functional enrichment disabled in config: step skipped")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
})
load_project_functions()

set.seed(cfg$project$seed)
processed_dir <- ensure_dir(cfg$paths$processed)
tables_dir <- ensure_dir(cfg$paths$tables)
figures_dir <- ensure_dir(cfg$paths$figures)
enrichment_cfg <- cfg$enrichment
keytype <- cfg$annotation$keytype
dpi <- cfg$plots$dpi

de_table <- readRDS(file.path(processed_dir, "de_results.rds"))
orgdb <- load_orgdb(cfg$annotation$orgdb)

as_keys <- function(ids) {
  unique(if (keytype == "ENSEMBL") strip_ensembl_version(ids) else ids)
}

# 1. ORA de genes sobreexpresados y subexpresados por separado ----
universe <- as_keys(de_table$gene_id[!is.na(de_table$padj)])
gene_sets <- list(
  up = as_keys(de_table$gene_id[de_table$status == "Up"]),
  down = as_keys(de_table$gene_id[de_table$status == "Down"])
)
direction_titles <- c(up = "GO ORA: genes sobreexpresados",
                      down = "GO ORA: genes subexpresados")

ora <- lapply(names(gene_sets), function(direction) {
  log_step("GO ORA (", enrichment_cfg$ontology, ") for ",
           length(gene_sets[[direction]]), " ", direction, "-regulated genes")
  run_go_ora(gene_sets[[direction]], universe, orgdb, keytype, enrichment_cfg)
})
names(ora) <- names(gene_sets)

for (direction in names(ora)) {
  write_table(enrichment_to_tibble(ora[[direction]]), tables_dir,
              sprintf("go_ora_%s.tsv", direction))
  save_plot(plot_enrichment_dotplot(ora[[direction]],
                                    direction_titles[[direction]]),
            figures_dir, sprintf("go_ora_%s_dotplot.png", direction),
            width = 9, height = 7, dpi = dpi)
}

# 2. GSEA ----
gene_list <- ranked_gene_list(de_table, keytype)
log_step("GO GSEA on ", length(gene_list), " ranked genes")
gsea <- run_go_gsea(gene_list, orgdb, keytype, enrichment_cfg)

write_table(enrichment_to_tibble(gsea), tables_dir, "go_gsea.tsv")
save_plot(plot_gsea_dotplot(gsea), figures_dir, "go_gsea_dotplot.png",
          width = 9, height = 9, dpi = dpi)

if (has_terms(gsea)) {
  save_base_plot(
    function() print(gseaplot2(gsea, geneSetID = seq_len(min(3, nrow(gsea))),
                               pvalue_table = FALSE)),
    figures_dir, "go_gsea_running_score.png", width = 8, height = 6, dpi = dpi
  )
}

saveRDS(list(ora = ora, gsea = gsea),
        file.path(processed_dir, "enrichment.rds"))
