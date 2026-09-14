# Funciones para preparar metadatos, filtrar y anotar genes y dar formato a los
# resultados de DESeq2.

#' Variables de una fórmula de diseño escrita como texto
design_variables <- function(formula_string) {
  all.vars(stats::as.formula(formula_string))
}

#' Prepara los metadatos de las muestras para DESeq2
#'
#' Convierte las variables del diseño en factores y fija el nivel de
#' referencia, para que el fold change sea tratamiento frente a referencia.
#'
#' @param samples Data frame (o DataFrame) con los metadatos.
#' @param design Sección `design` de la configuración.
prepare_sample_metadata <- function(samples, design) {
  vars <- design_variables(design$formula)
  missing <- setdiff(vars, colnames(samples))
  if (length(missing) > 0) {
    stop("Design variables not found in sample metadata: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!design$factor %in% vars) {
    stop("The contrast factor '", design$factor,
         "' must be part of the design formula.", call. = FALSE)
  }

  for (var in vars) {
    if (is.character(samples[[var]]) || is.logical(samples[[var]])) {
      samples[[var]] <- factor(samples[[var]])
    }
  }

  contrast_factor <- droplevels(factor(samples[[design$factor]]))
  absent <- setdiff(c(design$reference, design$treatment),
                    levels(contrast_factor))
  if (length(absent) > 0) {
    stop("Levels not found in '", design$factor, "': ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  samples[[design$factor]] <- stats::relevel(contrast_factor,
                                             ref = design$reference)
  samples
}

#' Indica qué genes tienen conteos suficientes para analizarse
#'
#' @param counts Matriz de conteos.
#' @param min_count Conteo mínimo por muestra.
#' @param min_samples Nº mínimo de muestras que alcanzan `min_count`.
#' @return Vector lógico con un valor por gen.
keep_expressed_genes <- function(counts, min_count = 10, min_samples = 3) {
  rowSums(counts >= min_count) >= min_samples
}

#' Aplica la transformación estabilizadora de la varianza (VST)
#'
#' `vst()` usa una submuestra de 1000 genes; con menos genes se aplica la
#' transformación completa.
variance_stabilize <- function(dds, blind = TRUE) {
  if (nrow(dds) >= 1000) {
    DESeq2::vst(dds, blind = blind)
  } else {
    DESeq2::varianceStabilizingTransformation(dds, blind = blind)
  }
}

#' Carga un paquete de anotación OrgDb por su nombre
load_orgdb <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Annotation package '", package, "' is not installed.",
         call. = FALSE)
  }
  getExportedValue(package, package)
}

#' Quita la versión de los IDs de Ensembl (ENSG00000000003.15 -> ENSG00000000003)
strip_ensembl_version <- function(ids) {
  sub("\\.[0-9]+$", "", ids)
}

#' Obtiene el símbolo, el nombre y el ID de Entrez de cada gen
#'
#' @param gene_ids Identificadores tal como aparecen en la matriz de conteos.
#' @param orgdb Objeto OrgDb, p. ej. org.Hs.eg.db.
#' @param keytype Tipo de identificador de `gene_ids`.
annotate_genes <- function(gene_ids, orgdb, keytype = "ENSEMBL") {
  keys <- if (keytype == "ENSEMBL") strip_ensembl_version(gene_ids) else gene_ids
  lookup <- function(column) {
    unname(suppressMessages(AnnotationDbi::mapIds(
      orgdb, keys = keys, keytype = keytype, column = column,
      multiVals = "first"
    )))
  }
  tibble::tibble(
    gene_id = gene_ids,
    symbol = lookup("SYMBOL"),
    gene_name = lookup("GENENAME"),
    entrez_id = lookup("ENTREZID")
  )
}

#' Clasifica los genes en sobreexpresados, subexpresados o no significativos
#'
#' @return Factor con los niveles Up, Down y NS.
classify_de <- function(log2fc, padj, alpha = 0.05, lfc_threshold = 1) {
  significant <- !is.na(padj) & !is.na(log2fc) &
    padj < alpha & abs(log2fc) >= lfc_threshold
  status <- ifelse(significant, ifelse(log2fc > 0, "Up", "Down"), "NS")
  factor(status, levels = c("Up", "Down", "NS"))
}

#' Une los resultados de DESeq2 en una tabla ordenada y anotada
#'
#' @param res `DESeqResults` de `results()` (fold change MLE, test de Wald).
#' @param res_shrunk `DESeqResults` de `lfcShrink()` (o el propio `res`).
#' @param annotation Tabla de `annotate_genes()`, o NULL.
#' @return Tibble ordenado por p-valor ajustado.
format_de_results <- function(res, res_shrunk, annotation = NULL,
                              alpha = 0.05, lfc_threshold = 1) {
  stopifnot(identical(rownames(res), rownames(res_shrunk)))

  de <- tibble::tibble(
    gene_id = rownames(res),
    baseMean = res$baseMean,
    log2FoldChange = res_shrunk$log2FoldChange,
    lfcSE = res_shrunk$lfcSE,
    log2FoldChange_mle = res$log2FoldChange,
    stat = res$stat,
    pvalue = res$pvalue,
    padj = res$padj
  )
  if (!is.null(annotation)) {
    de <- dplyr::left_join(de, annotation, by = "gene_id")
  }

  de |>
    dplyr::relocate(dplyr::any_of(c("symbol", "gene_name", "entrez_id")),
                    .after = "gene_id") |>
    dplyr::mutate(status = classify_de(.data$log2FoldChange, .data$padj,
                                       alpha, lfc_threshold)) |>
    dplyr::arrange(.data$padj, .data$pvalue)
}

#' Resume el número de genes diferencialmente expresados
summarise_de <- function(de_table) {
  tibble::tibble(
    genes_tested = nrow(de_table),
    genes_with_padj = sum(!is.na(de_table$padj)),
    up = sum(de_table$status == "Up"),
    down = sum(de_table$status == "Down"),
    total_significant = up + down
  )
}

#' Genes significativos con menor p-valor ajustado
top_de_genes <- function(de_table, n) {
  significant <- de_table[de_table$status != "NS", ]
  utils::head(significant[order(significant$padj), ], n)
}

#' Etiquetas de gen legibles: el símbolo si existe y, si no, el ID
gene_labels <- function(de_table) {
  if (!"symbol" %in% names(de_table)) {
    return(de_table$gene_id)
  }
  ifelse(is.na(de_table$symbol) | de_table$symbol == "",
         de_table$gene_id, de_table$symbol)
}
