# Enriquecimiento funcional con Gene Ontology: ORA y GSEA.

#' TRUE si el resultado de enriquecimiento tiene al menos un término
has_terms <- function(x) {
  !is.null(x) && nrow(as.data.frame(x)) > 0
}

#' Convierte un resultado de enriquecimiento en tibble (vacío si es NULL)
enrichment_to_tibble <- function(x) {
  if (is.null(x)) {
    return(tibble::tibble())
  }
  tibble::as_tibble(as.data.frame(x))
}

#' Análisis de sobrerrepresentación (ORA) de términos GO
#'
#' @param genes Identificadores de los genes de interés (p. ej. sobreexpresados).
#' @param universe Todos los genes analizados. Usarlos como fondo evita
#'   sobrestimar funciones que simplemente están expresadas.
#' @param enrichment_cfg Sección `enrichment` de la configuración.
#' @return Objeto `enrichResult`, o NULL si hay pocos genes.
run_go_ora <- function(genes, universe, orgdb, keytype, enrichment_cfg) {
  if (length(genes) < enrichment_cfg$min_genes) {
    log_step("Only ", length(genes), " genes provided: ORA skipped")
    return(NULL)
  }
  clusterProfiler::enrichGO(
    gene = genes,
    universe = universe,
    OrgDb = orgdb,
    keyType = keytype,
    ont = enrichment_cfg$ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = enrichment_cfg$pvalue_cutoff,
    qvalueCutoff = enrichment_cfg$qvalue_cutoff,
    minGSSize = enrichment_cfg$min_gs_size,
    maxGSSize = enrichment_cfg$max_gs_size,
    readable = TRUE
  )
}

#' Análisis GSEA de términos GO sobre una lista ordenada de genes
#'
#' @param gene_list Vector numérico ordenado de mayor a menor, con los IDs de
#'   gen como nombres.
#' @return Objeto `gseaResult` con símbolos de gen.
run_go_gsea <- function(gene_list, orgdb, keytype, enrichment_cfg) {
  gsea <- clusterProfiler::gseGO(
    geneList = gene_list,
    OrgDb = orgdb,
    keyType = keytype,
    ont = enrichment_cfg$ontology,
    minGSSize = enrichment_cfg$min_gs_size,
    maxGSSize = enrichment_cfg$max_gs_size,
    pvalueCutoff = enrichment_cfg$pvalue_cutoff,
    pAdjustMethod = "BH",
    eps = 0,
    seed = TRUE,
    verbose = FALSE
  )
  if (has_terms(gsea)) {
    gsea <- clusterProfiler::setReadable(gsea, OrgDb = orgdb,
                                         keyType = keytype)
  }
  gsea
}

#' Ordena los genes por el estadístico de Wald para el GSEA
#'
#' Este estadístico combina tamaño del efecto y precisión, y conserva el signo
#' del fold change.
ranked_gene_list <- function(de_table, keytype = "ENSEMBL") {
  ranked <- de_table[!is.na(de_table$stat), c("gene_id", "stat")]
  if (keytype == "ENSEMBL") {
    ranked$gene_id <- strip_ensembl_version(ranked$gene_id)
  }
  ranked <- ranked[!duplicated(ranked$gene_id), ]
  ranked <- ranked[order(ranked$stat, decreasing = TRUE), ]
  stats::setNames(ranked$stat, ranked$gene_id)
}
