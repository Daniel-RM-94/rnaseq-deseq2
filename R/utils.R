# Funciones generales: configuración, mensajes de registro, rutas y tablas.

#' Carga y valida la configuración del proyecto
#'
#' @param path Ruta al archivo YAML de configuración.
#' @return Lista con la configuración.
load_config <- function(path = here::here("config", "config.yml")) {
  if (!file.exists(path)) {
    stop("Configuration file not found: ", path, call. = FALSE)
  }
  cfg <- yaml::read_yaml(path)

  required <- c("project", "input", "design", "filtering", "de",
                "annotation", "enrichment", "plots", "paths")
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0) {
    stop("Missing configuration sections: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!cfg$input$source %in% c("airway", "files")) {
    stop("`input$source` must be either 'airway' or 'files'.", call. = FALSE)
  }
  if (!cfg$de$shrinkage %in% c("apeglm", "normal", "none")) {
    stop("`de$shrinkage` must be one of 'apeglm', 'normal' or 'none'.",
         call. = FALSE)
  }
  cfg
}

#' Carga todos los archivos de funciones de R/
load_project_functions <- function() {
  files <- list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)
  invisible(lapply(files, source))
}

#' Crea un directorio si no existe y devuelve su ruta absoluta
#'
#' @param path Ruta relativa a la raíz del proyecto.
ensure_dir <- function(path) {
  full_path <- here::here(path)
  dir.create(full_path, recursive = TRUE, showWarnings = FALSE)
  full_path
}

#' Muestra un mensaje con fecha y hora
log_step <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  paste0(...)))
}

#' Guarda un data frame como archivo separado por tabuladores
write_table <- function(x, dir, filename) {
  path <- file.path(dir, filename)
  readr::write_tsv(x, path, na = "")
  log_step("Saved table: ", path)
  invisible(path)
}

#' Lee una matriz de conteos (genes en filas, muestras en columnas)
#'
#' La primera columna debe tener identificadores de gen únicos.
#'
#' @param path Ruta a un archivo TSV, comprimido o no.
#' @return Matriz numérica de conteos.
read_count_matrix <- function(path) {
  df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  ids <- df[[1]]
  if (anyDuplicated(ids)) {
    stop("Duplicated gene identifiers found in ", path, call. = FALSE)
  }
  counts <- as.matrix(df[, -1, drop = FALSE])
  rownames(counts) <- ids
  counts
}

#' Lee la hoja de muestras (una fila por muestra)
#'
#' @param path Ruta a un archivo TSV.
#' @param id_column Columna con los identificadores de muestra, que deben
#'   coincidir con las columnas de la matriz de conteos.
read_sample_sheet <- function(path, id_column = "sample_id") {
  samples <- as.data.frame(
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  )
  if (!id_column %in% names(samples)) {
    stop("Column '", id_column, "' not found in ", path, call. = FALSE)
  }
  if (anyDuplicated(samples[[id_column]])) {
    stop("Duplicated sample identifiers found in ", path, call. = FALSE)
  }
  rownames(samples) <- samples[[id_column]]
  samples
}

#' Comprueba que los conteos y la hoja de muestras son coherentes
#'
#' @param counts Matriz de conteos crudos.
#' @param samples Data frame con los identificadores de muestra como nombres de fila.
validate_count_inputs <- function(counts, samples) {
  if (!is.numeric(counts)) {
    stop("The count matrix must be numeric.", call. = FALSE)
  }
  if (anyNA(counts)) {
    stop("The count matrix contains missing values.", call. = FALSE)
  }
  if (any(counts < 0)) {
    stop("The count matrix contains negative values.", call. = FALSE)
  }
  if (any(counts != round(counts))) {
    stop("The count matrix must contain raw integer counts ",
         "(not TPM/FPKM or other normalised values).", call. = FALSE)
  }
  if (!setequal(colnames(counts), rownames(samples))) {
    stop("Sample identifiers in the count matrix and the sample sheet ",
         "do not match.", call. = FALSE)
  }
  invisible(TRUE)
}
