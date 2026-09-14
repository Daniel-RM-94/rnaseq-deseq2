# Figuras de control de calidad, expresión diferencial y enriquecimiento.
# Las funciones con ggplot devuelven el gráfico; las de pheatmap lo dibujan o,
# si se indica `filename`, lo guardan en ese archivo.

DE_COLORS <- c(Up = "#B2182B", Down = "#2166AC", NS = "grey75")
DE_LABELS <- c("Sobreexpresados", "Subexpresados", "No significativos")

#' Tema de ggplot común a todo el proyecto
theme_rnaseq <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey30"),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      legend.key = ggplot2::element_blank()
    )
}

#' Guarda un ggplot en disco (no hace nada si el gráfico es NULL)
save_plot <- function(plot, dir, filename, width = 7, height = 5, dpi = 300) {
  if (is.null(plot)) {
    return(invisible(NULL))
  }
  path <- file.path(dir, filename)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi,
                  bg = "white")
  log_step("Saved figure: ", path)
  invisible(path)
}

#' Guarda en PNG un gráfico base o grid dibujado por `draw`
save_base_plot <- function(draw, dir, filename, width = 7, height = 5,
                           dpi = 300) {
  path <- file.path(dir, filename)
  grDevices::png(path, width = width, height = height, units = "in",
                 res = dpi)
  on.exit(grDevices::dev.off())
  draw()
  log_step("Saved figure: ", path)
  invisible(path)
}

# Control de calidad ----

#' Conteos totales por muestra
plot_library_sizes <- function(dds, fill_by) {
  df <- data.frame(
    sample = colnames(dds),
    millions = colSums(DESeq2::counts(dds)) / 1e6,
    group = SummarizedExperiment::colData(dds)[[fill_by]]
  )
  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(sample, millions),
                                   y = millions, fill = group)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_hline(yintercept = mean(df$millions), linetype = "dashed",
                        colour = "grey40") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(x = NULL, y = "Conteos totales (millones)", fill = fill_by,
                  title = "Tamaño de librería por muestra",
                  subtitle = "La línea discontinua indica la media") +
    theme_rnaseq()
}

#' Distribución de conteos crudos y normalizados por muestra
plot_count_distribution <- function(dds, fill_by) {
  to_long <- function(mat, label) {
    data.frame(
      sample = rep(colnames(mat), each = nrow(mat)),
      value = log2(as.vector(mat) + 1),
      type = label
    )
  }
  df <- rbind(
    to_long(DESeq2::counts(dds, normalized = FALSE), "Conteos crudos"),
    to_long(DESeq2::counts(dds, normalized = TRUE), "Conteos normalizados")
  )
  groups <- SummarizedExperiment::colData(dds)[[fill_by]]
  df$group <- groups[match(df$sample, colnames(dds))]

  ggplot2::ggplot(df, ggplot2::aes(x = sample, y = value, fill = group)) +
    ggplot2::geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
    ggplot2::facet_wrap(~type, ncol = 1) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(x = NULL, y = expression(log[2](conteos + 1)),
                  fill = fill_by,
                  title = "Distribución de conteos por muestra") +
    theme_rnaseq() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' PCA de los datos transformados con VST
plot_pca <- function(vsd, color_by, shape_by = NULL, ntop = 500) {
  pca <- DESeq2::plotPCA(vsd, intgroup = unique(c(color_by, shape_by)),
                         ntop = ntop, returnData = TRUE)
  percent_var <- 100 * attr(pca, "percentVar")

  mapping <- if (is.null(shape_by)) {
    ggplot2::aes(PC1, PC2, colour = .data[[color_by]])
  } else {
    ggplot2::aes(PC1, PC2, colour = .data[[color_by]],
                 shape = .data[[shape_by]])
  }

  ggplot2::ggplot(pca, mapping) +
    ggplot2::geom_point(size = 4, alpha = 0.9) +
    ggrepel::geom_text_repel(ggplot2::aes(label = name), size = 3,
                             colour = "grey30", show.legend = FALSE,
                             seed = 1) +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%% de la varianza)", percent_var[1]),
      y = sprintf("PC2 (%.1f%% de la varianza)", percent_var[2]),
      title = "Análisis de componentes principales",
      subtitle = sprintf("Datos VST, %d genes más variables", ntop)
    ) +
    theme_rnaseq()
}

#' Heatmap de distancias euclidianas entre muestras
plot_sample_distances <- function(vsd, annotation_cols, filename = NA) {
  distances <- stats::dist(t(SummarizedExperiment::assay(vsd)))
  annotation <- as.data.frame(
    SummarizedExperiment::colData(vsd)[, annotation_cols, drop = FALSE]
  )
  pheatmap::pheatmap(
    as.matrix(distances),
    clustering_distance_rows = distances,
    clustering_distance_cols = distances,
    annotation_col = annotation,
    color = grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(9, "Blues"))
    )(255),
    border_color = NA,
    main = "Distancia euclidiana entre muestras (VST)",
    filename = filename, width = 7, height = 6
  )
}

# Expresión diferencial ----

#' MA plot: expresión media frente a log2 fold change contraído
plot_ma <- function(de_table, lfc_threshold) {
  df <- de_table[order(de_table$status != "NS"), ]
  ggplot2::ggplot(df, ggplot2::aes(baseMean, log2FoldChange, colour = status)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey20") +
    ggplot2::geom_hline(yintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", colour = "grey40") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = DE_COLORS, breaks = names(DE_COLORS),
                                 labels = DE_LABELS, drop = FALSE) +
    ggplot2::labs(x = "Expresión media normalizada (baseMean)",
                  y = expression(log[2] ~ "fold change (contraído)"),
                  colour = NULL, title = "MA plot") +
    theme_rnaseq()
}

#' Volcano plot con los genes más significativos etiquetados
plot_volcano <- function(de_table, alpha, lfc_threshold, n_labels = 15) {
  df <- de_table[!is.na(de_table$padj), ]
  df$neg_log10_padj <- -log10(pmax(df$padj, .Machine$double.xmin))
  df <- df[order(df$status != "NS"), ]

  labels <- top_de_genes(df, n_labels)
  labels$label <- gene_labels(labels)

  counts <- c(sum(df$status == "Up"), sum(df$status == "Down"))
  legend_labels <- c(sprintf("%s (%d)", DE_LABELS[1:2], counts), DE_LABELS[3])

  ggplot2::ggplot(df, ggplot2::aes(log2FoldChange, neg_log10_padj,
                                   colour = status)) +
    ggplot2::geom_point(size = 1, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", colour = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed",
                        colour = "grey40") +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = label),
                             size = 3, colour = "black", max.overlaps = Inf,
                             min.segment.length = 0, seed = 1) +
    ggplot2::scale_colour_manual(values = DE_COLORS, breaks = names(DE_COLORS),
                                 labels = legend_labels, drop = FALSE) +
    ggplot2::labs(x = expression(log[2] ~ "fold change (contraído)"),
                  y = expression(-log[10] ~ "(padj)"), colour = NULL,
                  title = "Volcano plot",
                  subtitle = sprintf("Umbrales: padj < %s y |log2FC| ≥ %s",
                                     alpha, lfc_threshold)) +
    theme_rnaseq()
}

#' Heatmap de los genes DE más significativos (VST centrado por gen)
plot_top_genes_heatmap <- function(vsd, de_table, n = 50, annotation_cols,
                                   filename = NA) {
  top <- top_de_genes(de_table, n)
  if (nrow(top) < 2) {
    log_step("Fewer than 2 significant genes: heatmap skipped")
    return(invisible(NULL))
  }
  mat <- SummarizedExperiment::assay(vsd)[top$gene_id, , drop = FALSE]
  mat <- mat - rowMeans(mat)
  rownames(mat) <- make.unique(gene_labels(top))
  limit <- max(abs(mat))

  annotation <- as.data.frame(
    SummarizedExperiment::colData(vsd)[, annotation_cols, drop = FALSE]
  )
  pheatmap::pheatmap(
    mat,
    annotation_col = annotation,
    color = grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(11, "RdBu"))
    )(101),
    breaks = seq(-limit, limit, length.out = 102),
    border_color = NA,
    fontsize_row = 6,
    main = sprintf("Top %d genes DE (VST centrado por gen)", nrow(top)),
    filename = filename, width = 7, height = 9
  )
}

#' Conteos normalizados de los genes DE más significativos
plot_gene_counts <- function(dds, de_table, factor, color_by, n = 6) {
  top <- top_de_genes(de_table, n)
  if (nrow(top) == 0) {
    return(NULL)
  }
  top$label <- make.unique(gene_labels(top))

  df <- do.call(rbind, lapply(seq_len(nrow(top)), function(i) {
    d <- DESeq2::plotCounts(dds, gene = top$gene_id[i],
                            intgroup = unique(c(factor, color_by)),
                            returnData = TRUE)
    d$gene <- top$label[i]
    d
  }))
  df$gene <- base::factor(df$gene, levels = top$label)

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[factor]], y = count,
                                   colour = .data[[color_by]])) +
    ggplot2::geom_point(size = 2.5,
                        position = ggplot2::position_jitter(width = 0.15,
                                                            height = 0,
                                                            seed = 1)) +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~gene, scales = "free_y") +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::labs(x = NULL, y = "Conteos normalizados (escala log10)",
                  colour = color_by,
                  title = sprintf("Top %d genes diferencialmente expresados",
                                  nrow(top))) +
    theme_rnaseq()
}

# Enriquecimiento funcional ----

#' Dot plot de un resultado ORA
plot_enrichment_dotplot <- function(x, title, n = 15) {
  if (!has_terms(x)) {
    return(NULL)
  }
  enrichplot::dotplot(x, showCategory = n, label_format = 50) +
    ggplot2::ggtitle(title)
}

#' Dot plot del GSEA separado en términos activados y reprimidos
plot_gsea_dotplot <- function(x, n = 10) {
  if (!has_terms(x)) {
    return(NULL)
  }
  enrichplot::dotplot(x, showCategory = n, split = ".sign",
                      label_format = 50) +
    ggplot2::facet_grid(
      .sign ~ ., scales = "free_y", space = "free_y",
      labeller = ggplot2::as_labeller(c(activated = "Activados",
                                        suppressed = "Reprimidos"))
    ) +
    ggplot2::ggtitle("GSEA (Gene Ontology)")
}
