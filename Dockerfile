# syntax=docker/dockerfile:1
# Entorno reproducible para el análisis RNA-seq con DESeq2.
# Parte de la imagen oficial de Bioconductor (R, Bioconductor y RStudio Server).

# La imagen base se fija por digest porque la etiqueta se actualiza en origen.
# Para obtener el digest actual:
#   docker buildx imagetools inspect bioconductor/bioconductor_docker:RELEASE_3_23
ARG BIOC_VERSION=RELEASE_3_23
ARG BIOC_DIGEST=sha256:821dbf9ac119eac41f177531c7ca8fc7084c99c04eb218a1aa7f52cb63bad5d9
FROM bioconductor/bioconductor_docker:${BIOC_VERSION}@${BIOC_DIGEST}

LABEL org.opencontainers.image.title="rnaseq-deseq2" \
      org.opencontainers.image.description="Reproducible RNA-seq differential expression analysis with DESeq2 and Bioconductor" \
      org.opencontainers.image.licenses="MIT"

ENV PROJECT_HOME=/home/rstudio/project
WORKDIR ${PROJECT_HOME}

# Los paquetes se instalan antes de copiar el código para aprovechar la caché.
COPY docker/install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R && rm -f /tmp/install_packages.R

# RStudio Server se abre directamente en la carpeta del proyecto.
RUN echo "session-default-working-dir=${PROJECT_HOME}" >> /etc/rstudio/rsession.conf

# Script de arranque que exige una contraseña segura.
COPY --chmod=755 docker/start_rstudio.sh /usr/local/bin/start-rstudio

COPY --chown=rstudio:rstudio . ${PROJECT_HOME}

EXPOSE 8787

# Por defecto arranca RStudio Server. Los servicios pipeline y tests
# sustituyen este comando.
CMD ["/usr/local/bin/start-rstudio"]
