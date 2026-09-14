# Metodología

Este documento describe las decisiones metodológicas del análisis y su
justificación estadística. Los parámetros concretos se definen en
[`config/config.yml`](../config/config.yml).

## 1. Datos de entrada

El análisis parte de una **matriz de conteos crudos por gen**, es decir, el
número de lecturas asignadas a cada gen en cada muestra. DESeq2 modela
explícitamente la naturaleza discreta de los conteos y su relación
media-varianza, por lo que **no deben usarse valores normalizados** como TPM,
FPKM o CPM.

El dataset por defecto es `airway` (Himes et al. 2014, GSE52778). Contiene 8
muestras de 4 líneas celulares de músculo liso de vía aérea humana, cada una
tratada con dexametasona 1 µM durante 18 h (`trt`) o sin tratar (`untrt`).

## 2. Filtrado de genes poco expresados

Se conservan los genes con **≥ 10 conteos en al menos 4 muestras**, donde 4 es
el tamaño del grupo más pequeño. Así, un gen expresado exclusivamente en una
de las condiciones no se descarta. El filtrado:

- reduce el tamaño de los objetos y el tiempo de cálculo;
- elimina genes sin potencia estadística para detectar cambios;
- mejora la estimación de la tendencia de dispersión.

DESeq2 aplica además un **filtrado independiente** automático en `results()`,
que optimiza el número de descubrimientos para el FDR elegido.

## 3. Normalización

Las diferencias de profundidad de secuenciación y de composición entre
librerías se corrigen con los **factores de tamaño** de DESeq2 (método de la
*mediana de cocientes*). Para cada gen se calcula la media geométrica entre
muestras. El factor de cada muestra es la mediana de los cocientes entre sus
conteos y esas medias geométricas. El método es robusto frente a genes muy
expresados o diferencialmente expresados.

## 4. Modelo estadístico

Para cada gen *i* y muestra *j*, los conteos se modelan con una distribución
**binomial negativa**:

$$K_{ij} \sim \mathrm{NB}(\mu_{ij}, \alpha_i), \qquad \mu_{ij} = s_j \, q_{ij}, \qquad \log_2 q_{ij} = \sum_r x_{jr} \beta_{ir}$$

donde *s<sub>j</sub>* es el factor de tamaño, *α<sub>i</sub>* la dispersión del
gen y *β* los coeficientes del diseño.

### Diseño: `~ cell + dex`

Cada línea celular contribuye con una muestra tratada y otra control. Incluir
`cell` como covariable **bloquea** el efecto del donante: el coeficiente de
`dex` se estima a partir de las diferencias dentro de cada línea celular. Es el
equivalente, en un GLM, a un test pareado, y aumenta la potencia al eliminar la
variabilidad entre líneas de la varianza residual.

La variable de interés se sitúa al final de la fórmula por convención, y
`untrt` se fija como nivel de referencia. Por tanto, un log2FC positivo indica
mayor expresión con dexametasona.

### Estimación de la dispersión

Con pocas réplicas, la dispersión gen a gen es muy imprecisa. DESeq2 estima una
tendencia de la dispersión en función de la expresión media y **contrae** las
estimaciones individuales hacia esa tendencia mediante Bayes empírico. Los
genes con dispersión muy superior a la tendencia conservan su valor original
para no subestimar su variabilidad.

## 5. Test de hipótesis

- **Test de Wald** sobre el coeficiente `dex_trt_vs_untrt`.
- **Corrección por comparaciones múltiples** con Benjamini-Hochberg (FDR).
- Se usa `alpha = 0.05` tanto como umbral de FDR como objetivo del filtrado
  independiente.

## 6. Contracción del log2 fold change (apeglm)

Los fold changes de máxima verosimilitud (MLE) de genes con pocos conteos o
alta dispersión son muy ruidosos y pueden ser enormes sin ser fiables. El
método **apeglm** (Zhu et al. 2019) utiliza una distribución a priori de colas
pesadas (Cauchy) que:

- contrae fuertemente los fold changes poco informativos hacia cero;
- preserva los cambios grandes cuando la evidencia es sólida.

Los fold changes contraídos se usan para **ordenar, filtrar por magnitud y
visualizar** los resultados. Los p-valores provienen del test de Wald
original. Ambas versiones del fold change se incluyen en las tablas
(`log2FoldChange` y `log2FoldChange_mle`).

## 7. Criterio de gen diferencialmente expresado

Un gen se considera diferencialmente expresado (DE) si cumple:

- `padj < 0.05`, y
- `|log2FoldChange contraído| ≥ 1`, es decir, un cambio de al menos 2 veces.

Combinar significancia estadística y magnitud del efecto evita declarar
relevantes cambios diminutos que solo son significativos por su baja varianza.

> **Alternativa más estricta:** `results(dds, lfcThreshold = 1)` contrasta
> H0: |log2FC| ≤ 1 en lugar de H0: log2FC = 0. Esto controla formalmente el FDR
> sobre la magnitud del efecto, a costa de detectar menos genes.

## 8. Transformaciones para visualización

La PCA, los heatmaps y las distancias entre muestras usan la
**transformación estabilizadora de la varianza (VST)**. En la escala log
simple, la varianza de los genes poco expresados queda inflada; la VST corrige
ese efecto.

- **QC** (`blind = TRUE`): la transformación ignora el diseño, lo que permite
  una evaluación no sesgada de la estructura de los datos.
- **Visualización de genes DE** (`blind = FALSE`): usa la tendencia de
  dispersión del modelo para una transformación más precisa.

Los conteos transformados **solo se usan para explorar y visualizar**. El
análisis diferencial siempre se hace sobre los conteos crudos.

## 9. Anotación

Los identificadores Ensembl se traducen a símbolo, nombre y Entrez ID con
`AnnotationDbi::mapIds()` y `org.Hs.eg.db`. Si un identificador tiene varias
correspondencias, se toma la primera (`multiVals = "first"`).

## 10. Enriquecimiento funcional (Gene Ontology)

### ORA (*over-representation analysis*)

Se aplica un test hipergeométrico, con `clusterProfiler::enrichGO`, por
separado para los genes sobreexpresados y los subexpresados. Analizarlos por
separado evita que las señales de direcciones opuestas se diluyan.

- **Universo**: todos los genes con `padj` disponible, es decir, los genes que
  realmente podían detectarse como DE. Usar el genoma completo como fondo
  inflaría artificialmente el enriquecimiento de funciones generales de las
  células expresadas.
- Ontología: Proceso Biológico (BP); tamaño de conjuntos entre 10 y 500 genes;
  corrección BH.

### GSEA (*gene set enrichment analysis*)

Se usa `clusterProfiler::gseGO`, basado en `fgsea`, sobre **todos los genes**
ordenados por el **estadístico de Wald**. Este estadístico integra la magnitud
y la precisión del cambio y conserva su signo. A diferencia del ORA, el GSEA no
requiere umbrales arbitrarios y detecta cambios coordinados pero moderados en
conjuntos de genes.

- **NES > 0**: el conjunto de genes tiende a estar sobreexpresado con el
  tratamiento (activado).
- **NES < 0**: el conjunto tiende a estar subexpresado (reprimido).

## 11. Reproducibilidad

- Versiones fijadas mediante la imagen `bioconductor/bioconductor_docker:RELEASE_3_23`,
  referenciada por *digest* en el `Dockerfile` para que todos los builds partan
  exactamente de la misma imagen base. Una versión concreta de Bioconductor
  determina versiones compatibles de R y de todos los paquetes.
- Semilla aleatoria definida en la configuración (`project.seed`), utilizada
  por los pasos estocásticos (permutaciones de GSEA, submuestreo de `vst`).
- Cada paso del pipeline se ejecuta en una sesión de R independiente y se
  comunica con los demás solo a través de archivos. Así se evitan dependencias
  ocultas del entorno.
- `results/session_info.txt` y la última sección del informe registran las
  versiones exactas utilizadas.

## Referencias

- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550.
- Anders S, Huber W (2010). Differential expression analysis for sequence count data. *Genome Biology* 11:R106.
- Zhu A, Ibrahim JG, Love MI (2019). Heavy-tailed prior distributions for sequence count data. *Bioinformatics* 35(12):2084–2092.
- Himes BE et al. (2014). RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene. *PLoS ONE* 9(6):e99625.
- Wu T et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation* 2(3):100141.
- Korotkevich G et al. (2021). Fast gene set enrichment analysis. *bioRxiv* 060012.
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate. *J R Stat Soc B* 57(1):289–300.
- Bourgon R, Gentleman R, Huber W (2010). Independent filtering increases detection power for high-throughput experiments. *PNAS* 107(21):9546–9551.
