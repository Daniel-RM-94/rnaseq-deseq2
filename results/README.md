# Resultados

Esta carpeta se genera al ejecutar el pipeline y **no se versiona**, salvo este
README. En GitHub Actions los resultados se publican como *artifact*
(`rnaseq-results`) en cada ejecución.

```
results/
├── figures/
│   ├── qc_library_sizes.png          # Tamaño de librería por muestra
│   ├── qc_count_distribution.png     # Distribución de conteos crudos y normalizados
│   ├── qc_pca.png                    # PCA sobre datos VST
│   ├── qc_sample_distances.png       # Heatmap de distancias entre muestras
│   ├── de_dispersion.png             # Estimaciones de dispersión de DESeq2
│   ├── de_ma_plot.png                # MA plot
│   ├── de_volcano.png                # Volcano plot
│   ├── de_heatmap_top_genes.png      # Heatmap de los genes DE más significativos
│   ├── de_top_genes_counts.png       # Conteos normalizados de los genes top
│   ├── go_ora_{up,down}_dotplot.png  # GO ORA por dirección del cambio
│   ├── go_gsea_dotplot.png           # GO GSEA
│   └── go_gsea_running_score.png     # Curvas de enrichment score de GSEA
├── tables/
│   ├── sample_metadata.tsv
│   ├── filtering_summary.tsv
│   ├── qc_sample_metrics.tsv
│   ├── de_results_all.tsv            # Todos los genes analizados
│   ├── de_results_significant.tsv    # Solo genes DE
│   ├── de_summary.tsv
│   ├── go_ora_{up,down}.tsv
│   └── go_gsea.tsv
├── reports/
│   └── rnaseq_report.html            # Informe completo, autocontenido
└── session_info.txt                  # Versiones de R y de los paquetes
```

## Columnas de `de_results_*.tsv`

| Columna              | Descripción                                                     |
|----------------------|-----------------------------------------------------------------|
| `gene_id`            | Identificador Ensembl                                           |
| `symbol`             | Símbolo HGNC                                                    |
| `gene_name`          | Nombre completo del gen                                         |
| `entrez_id`          | Identificador Entrez                                            |
| `baseMean`           | Media de los conteos normalizados en todas las muestras         |
| `log2FoldChange`     | log2 fold change **contraído** (tratamiento vs referencia)      |
| `lfcSE`              | Error estándar del log2FC contraído                             |
| `log2FoldChange_mle` | log2 fold change de máxima verosimilitud, sin contraer          |
| `stat`               | Estadístico de Wald                                             |
| `pvalue`             | p-valor del test de Wald                                        |
| `padj`               | p-valor ajustado por Benjamini-Hochberg (vacío si el gen se excluyó en el filtrado independiente) |
| `status`             | `Up`, `Down` o `NS` según los umbrales de `config/config.yml`   |
