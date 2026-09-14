# Changelog

Todos los cambios relevantes se documentan en este archivo.
El formato sigue [Keep a Changelog](https://keepachangelog.com/es-ES/1.1.0/) y
el proyecto usa [versionado semántico](https://semver.org/lang/es/).

## [1.0.0] - 2026-09-13

### Añadido
- Pipeline modular de RNA-seq: importación, control de calidad, expresión
  diferencial con DESeq2 (contracción apeglm) y enriquecimiento GO (ORA y GSEA).
- Entorno reproducible con Docker (Bioconductor 3.23) y RStudio Server.
- Configuración centralizada en `config/config.yml` y soporte para datos propios.
- Informe HTML en R Markdown con interpretación de resultados.
- Tests unitarios con `testthat` e integración continua con GitHub Actions.

### Seguridad
- RStudio Server sin contraseña por defecto: `RSTUDIO_PASSWORD` debe definirse
  en `.env`. Se valida al arrancar (obligatoria, al menos 12 caracteres, no
  puede ser un valor por defecto y no se permite desactivar la autenticación).
- Puerto de RStudio Server enlazado solo a `127.0.0.1` y *healthcheck* del
  servicio.
- Política de seguridad en `SECURITY.md`.
