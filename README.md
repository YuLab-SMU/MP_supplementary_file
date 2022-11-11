<!-- README.md is generated from README.Rmd. Please edit that file -->

# MicrobiotaProcess: A comprehensive R package for managing and analyzing microbiome and other ecological data within the tidy framework

If you use this work in published research, please cite:

S Xu, L Zhan, W Tang, Z Dai, L Zhou, T Feng, M Chen, S Liu, X Fu, T Wu,
E Hu, G Yu. MicrobiotaProcess: A comprehensive R package for managing
and analyzing microbiome and other ecological data within the tidy
framework. 04 February 2022, [PREPRINT (Version 1) available at Research
Square](https://doi.org/10.21203/rs.3.rs-1284357/v1).

This repo contains source code and data to produce Supplementary
Material of the above paper.

To compile the
[supplementary\_file\_A.pdf](https://github.com/YuLab-SMU/MP_supplementary_file/blob/main/supplemental_file_A.pdf)
and
[supplementary\_file\_B.pdf](https://github.com/YuLab-SMU/MP_supplementary_file/blob/main/supplemental_file_B.pdf),
please run the following command in R:

    rmarkdown::render("supplementary_file_A.Rmd")
    rmarkdown::render("supplementary_file_B.Rmd")

Here is the output of `sessioninfo::session_info()` on the system on
which [this
document](https://github.com/YuLab-SMU/MP_supplementary_file/blob/main/supplemental_file.pdf)
was compiled:

    ## ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.2.0 (2022-04-22)
    ##  os       Ubuntu 18.04.4 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       Asia/Shanghai
    ##  date     2022-11-11
    ##  pandoc   2.9.2 @ /usr/bin/ (via rmarkdown)
    ## 
    ## ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  package                  * version    date (UTC) lib source
    ##  ade4                       1.7-19     2022-04-19 [1] CRAN (R 4.2.0)
    ##  ANCOMBC                  * 1.6.2      2022-06-26 [1] Bioconductor
    ##  AnnotationDbi              1.58.0     2022-04-26 [1] Bioconductor
    ##  AnnotationHub              3.4.0      2022-04-26 [1] Bioconductor
    ##  ape                        5.6-3      2022-10-30 [1] Github (emmanuelparadis/ape@090e82c)
    ##  aplot                    * 0.1.6.001  2022-08-26 [1] local
    ##  assertthat                 0.2.1      2019-03-21 [1] CRAN (R 4.2.0)
    ##  attempt                    0.3.1      2020-05-03 [1] CRAN (R 4.2.0)
    ##  backports                  1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
    ##  base64enc                  0.1-3      2015-07-28 [1] CRAN (R 4.2.0)
    ##  beachmat                   2.12.0     2022-04-26 [1] Bioconductor
    ##  beeswarm                   0.4.0      2021-06-01 [1] CRAN (R 4.2.0)
    ##  Biobase                  * 2.56.0     2022-04-26 [1] Bioconductor
    ##  BiocFileCache              2.4.0      2022-04-26 [1] Bioconductor
    ##  BiocGenerics             * 0.42.0     2022-04-26 [1] Bioconductor
    ##  BiocManager                1.30.18    2022-05-18 [1] CRAN (R 4.2.0)
    ##  BiocNeighbors              1.14.0     2022-04-26 [1] Bioconductor
    ##  BiocParallel               1.30.3     2022-06-05 [1] Bioconductor
    ##  BiocSingular               1.12.0     2022-04-26 [1] Bioconductor
    ##  BiocVersion                3.15.2     2022-03-29 [1] Bioconductor
    ##  biomformat                 1.24.0     2022-04-26 [1] Bioconductor
    ##  Biostrings               * 2.64.1     2022-08-18 [1] Bioconductor
    ##  bit                        4.0.4      2020-08-04 [1] CRAN (R 4.2.0)
    ##  bit64                      4.0.5      2020-08-30 [1] CRAN (R 4.2.0)
    ##  bitops                     1.0-7      2021-04-24 [1] CRAN (R 4.2.0)
    ##  blob                       1.2.3      2022-04-10 [1] CRAN (R 4.2.0)
    ##  boot                       1.3-28     2021-05-03 [1] CRAN (R 4.2.0)
    ##  broom                    * 1.0.0      2022-07-01 [1] CRAN (R 4.2.0)
    ##  bslib                      0.4.0      2022-07-16 [1] CRAN (R 4.2.0)
    ##  cachem                     1.0.6      2021-08-19 [1] CRAN (R 4.2.0)
    ##  caTools                    1.18.2     2021-03-28 [1] CRAN (R 4.2.0)
    ##  cellranger                 1.1.0      2016-07-27 [1] CRAN (R 4.2.0)
    ##  checkmate                  2.1.0      2022-04-21 [1] CRAN (R 4.2.0)
    ##  class                      7.3-20     2022-01-16 [1] CRAN (R 4.2.0)
    ##  cli                        3.4.1      2022-09-23 [1] CRAN (R 4.2.0)
    ##  clue                       0.3-61     2022-05-30 [1] CRAN (R 4.2.0)
    ##  cluster                    2.1.3      2022-03-28 [1] CRAN (R 4.2.0)
    ##  clusterProfiler          * 4.5.2      2022-09-06 [1] Bioconductor
    ##  codetools                  0.2-18     2020-11-04 [1] CRAN (R 4.2.0)
    ##  coin                     * 1.4-2      2021-10-08 [1] CRAN (R 4.2.0)
    ##  colorspace                 2.0-3      2022-02-21 [1] CRAN (R 4.2.0)
    ##  config                     0.3.1      2020-12-17 [1] CRAN (R 4.2.0)
    ##  crayon                     1.5.1      2022-03-26 [1] CRAN (R 4.2.0)
    ##  curatedMetagenomicData   * 3.4.2      2022-05-19 [1] Bioconductor
    ##  curl                       4.3.2      2021-06-23 [1] CRAN (R 4.2.0)
    ##  data.table                 1.14.2     2021-09-27 [1] CRAN (R 4.2.0)
    ##  DBI                        1.1.3      2022-06-18 [1] CRAN (R 4.2.0)
    ##  dbplyr                     2.2.1      2022-06-27 [1] CRAN (R 4.2.0)
    ##  DECIPHER                   2.24.0     2022-04-26 [1] Bioconductor
    ##  decontam                   1.16.0     2022-04-26 [1] Bioconductor
    ##  DelayedArray               0.22.0     2022-04-26 [1] Bioconductor
    ##  DelayedMatrixStats         1.18.0     2022-04-26 [1] Bioconductor
    ##  deldir                     1.0-6      2021-10-23 [1] CRAN (R 4.2.0)
    ##  desc                       1.4.2      2022-09-08 [1] CRAN (R 4.2.0)
    ##  DescTools                  0.99.46    2022-09-01 [1] CRAN (R 4.2.0)
    ##  digest                     0.6.30     2022-10-18 [1] CRAN (R 4.2.0)
    ##  DirichletMultinomial       1.38.0     2022-04-26 [1] Bioconductor
    ##  DO.db                      2.9        2022-09-06 [1] Bioconductor
    ##  doParallel                 1.0.17     2022-02-07 [1] CRAN (R 4.2.0)
    ##  doRNG                      1.8.2      2020-01-27 [1] CRAN (R 4.2.0)
    ##  DOSE                       3.23.2.001 2022-09-06 [1] Bioconductor
    ##  downloader                 0.4        2015-07-09 [1] CRAN (R 4.2.0)
    ##  dplyr                      1.0.10     2022-09-01 [1] CRAN (R 4.2.0)
    ##  DT                         0.25       2022-09-12 [1] CRAN (R 4.2.0)
    ##  e1071                      1.7-11     2022-06-07 [1] CRAN (R 4.2.0)
    ##  edgeR                    * 3.38.4     2022-08-07 [1] Bioconductor
    ##  ellipsis                   0.3.2      2021-04-29 [1] CRAN (R 4.2.0)
    ##  energy                     1.7-10     2022-04-19 [1] CRAN (R 4.2.0)
    ##  enrichplot               * 1.17.0.995 2022-09-06 [1] Bioconductor
    ##  evaluate                   0.16       2022-08-09 [1] CRAN (R 4.2.0)
    ##  Exact                      3.1        2021-11-26 [1] CRAN (R 4.2.0)
    ##  ExperimentHub              2.4.0      2022-04-26 [1] Bioconductor
    ##  expm                       0.999-6    2021-01-13 [1] CRAN (R 4.2.0)
    ##  fansi                      1.0.3      2022-03-24 [1] CRAN (R 4.2.0)
    ##  farver                     2.1.1      2022-07-06 [1] CRAN (R 4.2.0)
    ##  fastmap                    1.1.0      2021-01-25 [1] CRAN (R 4.2.0)
    ##  fastmatch                  1.1-3      2021-07-23 [1] CRAN (R 4.2.0)
    ##  fBasics                    4021.92    2022-08-08 [1] CRAN (R 4.2.0)
    ##  fgsea                      1.22.0     2022-04-26 [1] Bioconductor
    ##  filelock                   1.0.2      2018-10-05 [1] CRAN (R 4.2.0)
    ##  forcats                  * 0.5.1      2021-01-27 [1] CRAN (R 4.2.0)
    ##  foreach                    1.5.2      2022-02-02 [1] CRAN (R 4.2.0)
    ##  foreign                    0.8-82     2022-01-16 [1] CRAN (R 4.2.0)
    ##  Formula                    1.2-4      2020-10-16 [1] CRAN (R 4.2.0)
    ##  fs                         1.5.2      2021-12-08 [1] CRAN (R 4.2.0)
    ##  generics                   0.1.3      2022-07-05 [1] CRAN (R 4.2.0)
    ##  GenomeInfoDb             * 1.32.4     2022-09-06 [1] Bioconductor
    ##  GenomeInfoDbData           1.2.8      2022-04-28 [1] Bioconductor
    ##  GenomicRanges            * 1.48.0     2022-04-26 [1] Bioconductor
    ##  ggbeeswarm                 0.6.0      2017-08-07 [1] CRAN (R 4.2.0)
    ##  ggforce                    0.3.3      2021-03-05 [1] CRAN (R 4.2.0)
    ##  ggfun                      0.0.6      2022-08-30 [1] local
    ##  ggnewscale               * 0.4.7      2022-03-25 [1] CRAN (R 4.2.0)
    ##  ggplot2                  * 3.4.0      2022-11-04 [1] CRAN (R 4.2.0)
    ##  ggplotify                  0.1.0      2021-09-02 [1] CRAN (R 4.2.0)
    ##  ggpp                     * 0.4.5      2022-09-30 [1] CRAN (R 4.2.0)
    ##  ggraph                     2.0.6      2022-08-08 [1] CRAN (R 4.2.0)
    ##  ggrepel                  * 0.9.1      2021-01-15 [1] CRAN (R 4.2.0)
    ##  ggsignif                   0.6.3      2021-09-09 [1] CRAN (R 4.2.0)
    ##  ggstar                     1.0.3      2021-12-03 [1] CRAN (R 4.2.0)
    ##  ggtree                   * 3.7.1      2022-11-10 [1] Bioconductor
    ##  ggtreeExtra              * 1.8.0      2022-11-04 [1] Bioconductor
    ##  ggupset                  * 0.3.0      2020-05-05 [1] CRAN (R 4.2.0)
    ##  ggVennDiagram            * 1.2.0      2021-10-22 [1] CRAN (R 4.2.0)
    ##  gld                        2.6.5      2022-06-29 [1] CRAN (R 4.2.0)
    ##  glmnet                   * 4.1-4      2022-04-15 [1] CRAN (R 4.2.0)
    ##  glue                       1.6.2      2022-02-24 [1] CRAN (R 4.2.0)
    ##  GO.db                      3.15.0     2022-09-06 [1] Bioconductor
    ##  golem                      0.3.5      2022-10-18 [1] CRAN (R 4.2.0)
    ##  GOSemSim                   2.22.0     2022-04-26 [1] Bioconductor
    ##  gplots                     3.1.3      2022-04-25 [1] CRAN (R 4.2.0)
    ##  graphlayouts               0.8.0      2022-01-03 [1] CRAN (R 4.2.0)
    ##  gridExtra                  2.3        2017-09-09 [1] CRAN (R 4.2.0)
    ##  gridGraphics               0.5-1      2020-12-13 [1] CRAN (R 4.2.0)
    ##  gsl                        2.1-7.1    2021-11-02 [1] CRAN (R 4.2.0)
    ##  gson                       0.0.8      2022-08-20 [1] CRAN (R 4.2.0)
    ##  gtable                     0.3.1      2022-09-01 [1] CRAN (R 4.2.0)
    ##  gtools                     3.9.3      2022-07-11 [1] CRAN (R 4.2.0)
    ##  GUniFrac                 * 1.7        2022-10-23 [1] CRAN (R 4.2.0)
    ##  Hmisc                      4.7-0      2022-04-19 [1] CRAN (R 4.2.0)
    ##  hms                        1.1.1      2021-09-26 [1] CRAN (R 4.2.0)
    ##  htmlTable                  2.4.1      2022-07-07 [1] CRAN (R 4.2.0)
    ##  htmltools                  0.5.3      2022-07-18 [1] CRAN (R 4.2.0)
    ##  htmlwidgets                1.5.4      2021-09-08 [1] CRAN (R 4.2.0)
    ##  httpuv                     1.6.6      2022-09-08 [1] CRAN (R 4.2.0)
    ##  httr                       1.4.4      2022-08-17 [1] CRAN (R 4.2.0)
    ##  igraph                   * 1.3.4      2022-07-19 [1] CRAN (R 4.2.0)
    ##  interactiveDisplayBase     1.34.0     2022-04-26 [1] Bioconductor
    ##  interp                     1.1-3      2022-07-13 [1] CRAN (R 4.2.0)
    ##  IRanges                  * 2.30.1     2022-08-18 [1] Bioconductor
    ##  irlba                      2.3.5      2021-12-06 [1] CRAN (R 4.2.0)
    ##  iterators                  1.0.14     2022-02-05 [1] CRAN (R 4.2.0)
    ##  jpeg                       0.1-9      2021-07-24 [1] CRAN (R 4.2.0)
    ##  jquerylib                  0.1.4      2021-04-26 [1] CRAN (R 4.2.0)
    ##  jsonlite                   1.8.0      2022-02-22 [1] CRAN (R 4.2.0)
    ##  KEGGREST                   1.36.3     2022-07-12 [1] Bioconductor
    ##  KernSmooth                 2.23-20    2021-05-03 [1] CRAN (R 4.2.0)
    ##  knitr                      1.39       2022-04-26 [1] CRAN (R 4.2.0)
    ##  later                      1.3.0      2021-08-18 [1] CRAN (R 4.2.0)
    ##  lattice                    0.20-45    2021-09-22 [1] CRAN (R 4.2.0)
    ##  latticeExtra               0.6-30     2022-07-04 [1] CRAN (R 4.2.0)
    ##  lazyeval                   0.2.2      2019-03-15 [1] CRAN (R 4.2.0)
    ##  libcoin                    1.0-9      2021-09-27 [1] CRAN (R 4.2.0)
    ##  lifecycle                  1.0.3      2022-10-07 [1] CRAN (R 4.2.0)
    ##  limma                    * 3.52.3     2022-09-11 [1] Bioconductor
    ##  lme4                       1.1-30     2022-07-08 [1] CRAN (R 4.2.0)
    ##  lmerTest                   3.1-3      2020-10-23 [1] CRAN (R 4.2.0)
    ##  lmom                       2.9        2022-05-29 [1] CRAN (R 4.2.0)
    ##  locfit                     1.5-9.6    2022-07-11 [1] CRAN (R 4.2.0)
    ##  magrittr                   2.0.3      2022-03-30 [1] CRAN (R 4.2.0)
    ##  MASS                       7.3-57     2022-04-22 [1] CRAN (R 4.2.0)
    ##  Matrix                   * 1.4-1      2022-03-23 [1] CRAN (R 4.2.0)
    ##  MatrixGenerics           * 1.8.1      2022-06-26 [1] Bioconductor
    ##  matrixStats              * 0.62.0     2022-04-19 [1] CRAN (R 4.2.0)
    ##  memoise                    2.0.1      2021-11-26 [1] CRAN (R 4.2.0)
    ##  metagenomeSeq            * 1.38.0     2022-04-26 [1] Bioconductor
    ##  mgcv                       1.8-40     2022-03-29 [1] CRAN (R 4.2.0)
    ##  mia                        1.4.0      2022-04-26 [1] Bioconductor
    ##  microbiome                 1.18.0     2022-04-26 [1] Bioconductor
    ##  MicrobiomeProfiler       * 1.1.0      2022-11-08 [1] Bioconductor
    ##  MicrobiomeStat           * 1.1        2022-01-24 [1] CRAN (R 4.2.0)
    ##  MicrobiotaProcess        * 1.11.2     2022-11-09 [1] Bioconductor
    ##  mime                       0.12       2021-09-28 [1] CRAN (R 4.2.0)
    ##  minqa                      1.2.4      2014-10-09 [1] CRAN (R 4.2.0)
    ##  modeest                    2.4.0      2019-11-18 [1] CRAN (R 4.2.0)
    ##  modeltools                 0.2-23     2020-03-05 [1] CRAN (R 4.2.0)
    ##  multcomp                   1.4-19     2022-04-26 [1] CRAN (R 4.2.0)
    ##  MultiAssayExperiment       1.22.0     2022-04-26 [1] Bioconductor
    ##  multtest                   2.52.0     2022-04-26 [1] Bioconductor
    ##  munsell                    0.5.0      2018-06-12 [1] CRAN (R 4.2.0)
    ##  mvtnorm                    1.1-3      2021-10-08 [1] CRAN (R 4.2.0)
    ##  nlme                       3.1-158    2022-06-15 [1] CRAN (R 4.2.0)
    ##  nloptr                     2.0.3      2022-05-26 [1] CRAN (R 4.2.0)
    ##  nnet                       7.3-17     2022-01-16 [1] CRAN (R 4.2.0)
    ##  numDeriv                   2016.8-1.1 2019-06-06 [1] CRAN (R 4.2.0)
    ##  patchwork                * 1.1.2      2022-08-19 [1] CRAN (R 4.2.0)
    ##  permute                    0.9-7      2022-01-27 [1] CRAN (R 4.2.0)
    ##  phyloseq                   1.40.0     2022-04-26 [1] Bioconductor
    ##  pillar                     1.8.1      2022-08-19 [1] CRAN (R 4.2.0)
    ##  pkgconfig                  2.0.3      2019-09-22 [1] CRAN (R 4.2.0)
    ##  pkgload                    1.3.0      2022-06-27 [1] CRAN (R 4.2.0)
    ##  plyr                       1.8.7      2022-03-24 [1] CRAN (R 4.2.0)
    ##  png                        0.1-7      2013-12-03 [1] CRAN (R 4.2.0)
    ##  polyclip                   1.10-0     2019-03-14 [1] CRAN (R 4.2.0)
    ##  preprocessCore             1.58.0     2022-04-26 [1] Bioconductor
    ##  pROC                     * 1.18.0     2021-09-03 [1] CRAN (R 4.2.0)
    ##  promises                   1.2.0.1    2021-02-11 [1] CRAN (R 4.2.0)
    ##  proxy                      0.4-27     2022-06-09 [1] CRAN (R 4.2.0)
    ##  purrr                      0.3.5      2022-10-06 [1] CRAN (R 4.2.0)
    ##  qvalue                     2.28.0     2022-04-26 [1] Bioconductor
    ##  R6                         2.5.1      2021-08-19 [1] CRAN (R 4.2.0)
    ##  randomForest             * 4.7-1.1    2022-05-23 [1] CRAN (R 4.2.0)
    ##  rappdirs                   0.3.3      2021-01-31 [1] CRAN (R 4.2.0)
    ##  rbibutils                  2.2.9      2022-08-15 [1] CRAN (R 4.2.0)
    ##  RColorBrewer             * 1.1-3      2022-04-03 [1] CRAN (R 4.2.0)
    ##  Rcpp                       1.0.9      2022-07-08 [1] CRAN (R 4.2.0)
    ##  RCurl                      1.98-1.8   2022-07-30 [1] CRAN (R 4.2.0)
    ##  Rdpack                     2.4        2022-07-20 [1] CRAN (R 4.2.0)
    ##  readr                      2.1.2      2022-01-30 [1] CRAN (R 4.2.0)
    ##  readxl                     1.4.1      2022-08-17 [1] CRAN (R 4.2.0)
    ##  reshape2                   1.4.4      2020-04-09 [1] CRAN (R 4.2.0)
    ##  rhdf5                      2.40.0     2022-04-26 [1] Bioconductor
    ##  rhdf5filters               1.8.0      2022-04-26 [1] Bioconductor
    ##  Rhdf5lib                   1.18.2     2022-05-15 [1] Bioconductor
    ##  rlang                      1.0.6      2022-09-24 [1] CRAN (R 4.2.0)
    ##  rmarkdown                  2.15       2022-08-16 [1] CRAN (R 4.2.0)
    ##  rmutil                     1.1.9      2022-03-01 [1] CRAN (R 4.2.0)
    ##  rngtools                   1.5.2      2021-09-20 [1] CRAN (R 4.2.0)
    ##  rootSolve                  1.8.2.3    2021-09-29 [1] CRAN (R 4.2.0)
    ##  roxygen2                   7.2.1      2022-07-18 [1] CRAN (R 4.2.0)
    ##  rpart                      4.1.16     2022-01-24 [1] CRAN (R 4.2.0)
    ##  rprojroot                  2.0.3      2022-04-02 [1] CRAN (R 4.2.0)
    ##  RSQLite                    2.2.17     2022-09-10 [1] CRAN (R 4.2.0)
    ##  rstudioapi                 0.13       2020-11-12 [1] CRAN (R 4.2.0)
    ##  rsvd                       1.0.5      2021-04-16 [1] CRAN (R 4.2.0)
    ##  Rtsne                      0.16       2022-04-17 [1] CRAN (R 4.2.0)
    ##  RVenn                      1.1.0      2019-07-18 [1] CRAN (R 4.2.0)
    ##  S4Vectors                * 0.34.0     2022-04-26 [1] Bioconductor
    ##  sandwich                   3.0-2      2022-06-15 [1] CRAN (R 4.2.0)
    ##  sass                       0.4.2      2022-07-16 [1] CRAN (R 4.2.0)
    ##  ScaledMatrix               1.4.0      2022-04-26 [1] Bioconductor
    ##  scales                     1.2.1      2022-08-20 [1] CRAN (R 4.2.0)
    ##  scater                     1.24.0     2022-04-26 [1] Bioconductor
    ##  scatterpie                 0.1.7      2022-09-02 [1] local
    ##  scuttle                    1.6.3      2022-08-23 [1] Bioconductor
    ##  sessioninfo                1.2.2      2021-12-06 [1] CRAN (R 4.2.0)
    ##  shadowtext               * 0.1.2      2022-04-22 [1] CRAN (R 4.2.0)
    ##  shape                      1.4.6      2021-05-19 [1] CRAN (R 4.2.0)
    ##  shiny                      1.7.2      2022-07-19 [1] CRAN (R 4.2.0)
    ##  shinycustomloader          0.9.0      2018-03-27 [1] CRAN (R 4.2.0)
    ##  shinyWidgets               0.7.4      2022-10-05 [1] CRAN (R 4.2.0)
    ##  SingleCellExperiment     * 1.18.0     2022-04-26 [1] Bioconductor
    ##  sparseMatrixStats          1.8.0      2022-04-26 [1] Bioconductor
    ##  spatial                    7.3-15     2022-01-16 [1] CRAN (R 4.2.0)
    ##  stable                     1.1.6      2022-03-02 [1] CRAN (R 4.2.0)
    ##  stabledist                 0.7-1      2016-09-12 [1] CRAN (R 4.2.0)
    ##  statip                     0.2.3      2019-11-17 [1] CRAN (R 4.2.0)
    ##  statmod                    1.4.37     2022-08-12 [1] CRAN (R 4.2.0)
    ##  stringi                    1.7.8      2022-07-11 [1] CRAN (R 4.2.0)
    ##  stringr                    1.4.1      2022-08-20 [1] CRAN (R 4.2.0)
    ##  SummarizedExperiment     * 1.26.1     2022-04-29 [1] Bioconductor
    ##  survival                 * 3.3-1      2022-03-03 [1] CRAN (R 4.2.0)
    ##  TH.data                    1.1-1      2022-04-26 [1] CRAN (R 4.2.0)
    ##  tibble                     3.1.8      2022-07-22 [1] CRAN (R 4.2.0)
    ##  tidybulk                 * 1.8.2      2022-09-29 [1] Bioconductor
    ##  tidygraph                  1.2.2      2022-08-22 [1] CRAN (R 4.2.0)
    ##  tidyr                      1.2.1      2022-09-08 [1] CRAN (R 4.2.0)
    ##  tidyselect                 1.2.0      2022-10-10 [1] CRAN (R 4.2.0)
    ##  tidytree                   0.3.9      2022-09-19 [1] local
    ##  timeDate                   4021.104   2022-07-19 [1] CRAN (R 4.2.0)
    ##  timeSeries                 4021.104   2022-07-17 [1] CRAN (R 4.2.0)
    ##  treeio                     1.21.3     2022-10-30 [1] Bioconductor
    ##  TreeSummarizedExperiment * 2.4.0      2022-04-26 [1] Bioconductor
    ##  tweenr                     1.0.2      2021-03-23 [1] CRAN (R 4.2.0)
    ##  tzdb                       0.3.0      2022-03-28 [1] CRAN (R 4.2.0)
    ##  usethis                    2.1.6      2022-05-25 [1] CRAN (R 4.2.0)
    ##  utf8                       1.2.2      2021-07-24 [1] CRAN (R 4.2.0)
    ##  vctrs                      0.5.0      2022-10-22 [1] CRAN (R 4.2.0)
    ##  vegan                      2.6-2      2022-04-17 [1] CRAN (R 4.2.0)
    ##  vipor                      0.4.5      2017-03-22 [1] CRAN (R 4.2.0)
    ##  viridis                    0.6.2      2021-10-13 [1] CRAN (R 4.2.0)
    ##  viridisLite                0.4.1      2022-08-22 [1] CRAN (R 4.2.0)
    ##  withr                      2.5.0      2022-03-03 [1] CRAN (R 4.2.0)
    ##  Wrench                     1.14.0     2022-04-26 [1] Bioconductor
    ##  xfun                       0.32       2022-08-10 [1] CRAN (R 4.2.0)
    ##  xml2                       1.3.3      2021-11-30 [1] CRAN (R 4.2.0)
    ##  xtable                     1.8-4      2019-04-21 [1] CRAN (R 4.2.0)
    ##  XVector                  * 0.36.0     2022-04-26 [1] Bioconductor
    ##  yaml                       2.3.5      2022-02-21 [1] CRAN (R 4.2.0)
    ##  yulab.utils                0.0.5      2022-06-30 [1] CRAN (R 4.2.0)
    ##  zlibbioc                   1.42.0     2022-04-26 [1] Bioconductor
    ##  zoo                        1.8-10     2022-04-15 [1] CRAN (R 4.2.0)
    ## 
    ##  [1] /mnt/d/UbuntuApps/R/4.2.0/lib/R/library
    ## 
    ## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
