# Datos de entrada

Los archivos de esta carpeta **no se versionan** (ver `.gitignore`).

Con la configuración por defecto (`input.source: airway`), el pipeline carga el
dataset de ejemplo desde el paquete de Bioconductor `airway` y exporta aquí:

- `airway_counts.tsv.gz`: matriz de conteos.
- `airway_samples.tsv`: hoja de muestras.

Estos archivos sirven como referencia del formato esperado para tus propios datos.

## Usar tus propios datos

1. Coloca aquí dos archivos TSV (tabulados, opcionalmente comprimidos con gzip).

   **Matriz de conteos** (`counts.tsv.gz`): la primera columna contiene
   identificadores de gen únicos y cada columna siguiente corresponde a una
   muestra. Los valores deben ser **conteos crudos enteros**, por ejemplo de
   `featureCounts`, `HTSeq-count` o `tximport` con `countsFromAbundance`. No se
   admiten TPM ni FPKM.

   | gene_id         | sample_1 | sample_2 | sample_3 | sample_4 |
   |-----------------|---------:|---------:|---------:|---------:|
   | ENSG00000000003 |      679 |      448 |      873 |      408 |
   | ENSG00000000419 |      467 |      515 |      621 |      365 |

   **Hoja de muestras** (`samples.tsv`): una fila por muestra. La columna
   `sample_id` debe coincidir con las columnas de la matriz y debe haber una
   columna por cada variable del diseño.

   | sample_id | condition | batch |
   |-----------|-----------|-------|
   | sample_1  | control   | A     |
   | sample_2  | treated   | A     |

2. Edita `config/config.yml`:

   ```yaml
   input:
     source: files
     counts_file: data/raw/counts.tsv.gz
     samples_file: data/raw/samples.tsv
   design:
     formula: "~ batch + condition"
     factor: condition
     reference: control
     treatment: treated
   ```

3. Si tu organismo no es humano, cambia `annotation.orgdb` (por ejemplo
   `org.Mm.eg.db` para ratón), añade el paquete a `docker/install_packages.R` y
   reconstruye la imagen.
