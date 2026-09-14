# Instala los paquetes del proyecto. La imagen base fija la versión de
# Bioconductor, así que las versiones instaladas son coherentes entre sí.

packages <- c(
  # Bioconductor: expresión diferencial, anotación y enriquecimiento
  "DESeq2",
  "apeglm",
  "airway",
  "SummarizedExperiment",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichplot",
  # CRAN: manejo de datos, gráficos e informes
  "dplyr",
  "tibble",
  "readr",
  "yaml",
  "here",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "RColorBrewer",
  "rmarkdown",
  "knitr",
  "sessioninfo",
  "testthat"
)

BiocManager::install(
  packages,
  update = FALSE,
  ask = FALSE,
  Ncpus = max(1L, parallel::detectCores() - 1L)
)

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("The following packages failed to install: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

message("All ", length(packages), " packages installed (Bioconductor ",
        BiocManager::version(), ", ", R.version.string, ").")
