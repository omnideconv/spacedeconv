# :rocket: spacedeconv <a href="https://omnideconv.github.io/spacedeconv"><img src="man/figures/logo.png" align="right" height="200" /></a>

[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![R](https://img.shields.io/badge/R-4.3.3-276DC3?logo=r&logoColor=white)
![Python](https://img.shields.io/badge/Python-3.10-3776AB?logo=python&logoColor=white)
![Platform](https://img.shields.io/badge/platform-Linux%2064--bit-lightgrey)

_spacedeconv_ is a unified interface to first- and second-generation deconvolution tools with focus on spatial transcriptomics datasets. The package is able to directly estimate celltype proportions of immune cells and can deconvolute any celltype if an annotated single-cell reference dataset is available.

## :arrow_down: Installation

**Note:** The current _spacedeconv_ installation is only available for the **linux-64** platform.

Since many different packages need to be included, we highly recommend to install _spacedeconv_ in a new Conda environment with the following commands.

First, a tool for fast dependency resolution is needed, therefore we recommend installing mamba if not already available:

```r
conda install -c conda-forge mamba
```

Download the environment.yml file of this github repo:

```r
wget https://raw.githubusercontent.com/omnideconv/spacedeconv/main/environment.yml -O environment.yml
```

Create a new environment called "r-omnideconv" via mamba with the environment.yml file:

```r
mamba env create -f environment.yml
```

Start R inside the r-omnideconv conda environment:

```r
conda activate r-omnideconv
R
```

Install the missing packages that are not available via conda as well as _spacedeconv_:

```r
pak::pkg_install("drieslab/Giotto@v3.3.2", upgrade = FALSE)
devtools::install_github("YingMa0107/CARD", ref = "2d64b91abb5cdd0c7f576b1c5d4727c84e7c93a0", upgrade = "never")
pak::pkg_install("omnideconv/spacedeconv", dependencies = FALSE, upgrade = FALSE)
```

## :sparkles: Features

- unified access to first- and second-generation deconvolution tools
- direct deconvolution of immune cells
- compute custom reference signatures to deconvolute any celltype
- flexible visualization functions
- resource optimization
- Pathway and Transcription Factor analysis (decoupleR integration)
- Ligand-Receptor analysis
- easy integration into spatial transcriptomics workflows

## :floppy_disk: Data requirements

Spatial transcriptomics data: _[SpatialExperiment](https://bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html)_

Single-cell data with cell-type annotation: _[SingleCellExperiment](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)_ (recommended) or _[Seurat](https://satijalab.org/seurat/)_

## :technologist: Usage

The main workflow consists of the following steps:

### 1. Load datasets

To explore the package, start by loading some of the built-in example datasets.

```r
library(spacedeconv)

# data("single_cell_data_1")
# data("single_cell_data_2")
data("single_cell_data_3")
# data("single_cell_data_4")

# data("spatial_data_1")
# data("spatial_data_2")
data("spatial_data_3")
# data("spatial_data_4")
```

### 2. Preprocessing

Depending on the use case, certain preprocessing steps might be necessary.

```r
single_cell_data_3 <- spacedeconv::preprocess(single_cell_data_3)
spatial_data_3 <- spacedeconv::preprocess(spatial_data_3)

spatial_data_3 <- spacedeconv::normalize(spatial_data_3, method = "cpm")
```

### 3. Build a Signature Matrix

Build a cell type specific signature matrix from annotated single-cell reference data.

```r
signature <- spacedeconv::build_model(
  single_cell_obj = single_cell_data_3,
  cell_type_col = "celltype_major",
  method = "spatialdwls", verbose = T
)
```

### 4. Deconvolution

While some methods are able to directly estimate immune cell abundances other tools require a custom reference signature.

```r
deconv <- spacedeconv::deconvolute(
  spatial_obj = spatial_data_3,
  single_cell_obj = single_cell_data_3,
  cell_type_col = "celltype_major",
  method = "spatialdwls",
  signature = signature,
  assay_sp = "cpm"
)
```

### 5. Visualization

_spacedeconv_ includes highly-flexible visualization functions.

```r
plot_spatial(
  spe = deconv,
  result = "spatialdwls_B.cells",
  title = "B cells",
  density=F
)
```

## Available methods, Licenses, Citations

Note that, while _spacedeconv_ itself is free ([GPL
3.0](https://github.com/omnideconv/omnideconv/blob/main/LICENSE)), you may
need to obtain a license to use the individual methods. See the table
below for more information. If you use this package in your work, please
cite both our package and the method(s) you are using.

> Constantin Zackl, Maria Zopoglou, Reto Stauffer, Markus Ausserhofer, Marieke E. Ijsselsteijn, Gregor Sturm, Noel Filipe da Cunha Carvalho de Miranda, Francesca Finotello. spacedeconv: deconvolution of tissue architecture from spatial transcriptomics, PREPRINT available at Research Square https://doi.org/10.21203/rs.3.rs-5102166/v1.

| Method                                                         | Usable with signature() |                                     Licence                                     | Citation                                                                                                                                                                                                                            |
| -------------------------------------------------------------- | :---------------------: | :-----------------------------------------------------------------------------: | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [immunedeconv](https://github.com/omnideconv/immunedeconv)     |           :x:           |      [BSD](https://github.com/omnideconv/immunedeconv/blob/master/LICENSE)      | Sturm, G. et al. Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. Bioinformatics, 35(14), i436-i445 (2019). <https://doi.org/10.1093/bioinformatics/btz363>                    |
| [spatialDWLS](https://github.com/RubD/Giotto/)                 |   :heavy_check_mark:    |            [MIT](https://github.com/RubD/Giotto/blob/master/LICENSE)            | Dong R., Yuan G.C. SpatialDWLS: accurate deconvolution of spatial transcriptomic data. Genome Biology 22, 145 (2021). <https://doi.org/10.1186/s13059-021-02362-7>                                                                  |
| [cell2location](https://github.com/BayraktarLab/cell2location) |   :heavy_check_mark:    | [Apache-2.0](https://github.com/BayraktarLab/cell2location/blob/master/LICENSE) | Kleshchevnikov V. et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nature Biotechnology 40, 661–671 (2022). <https://doi.org/10.1038/s41587-021-01139-4>                                               |
| [SPOTlight](https://github.com/MarcElosua/SPOTlight)           |   :heavy_check_mark:    |     [GPL 3.0](https://github.com/MarcElosua/SPOTlight/blob/main/LICENSE.md)     | Elosua-Bayes M. et al. SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes. Nucleic Acids Research 49(9):e50 (2021). <https://doi.org/10.1093/nar/gkab043>.               |
| [RCTD](https://github.com/dmcable/spacexr)                     |   :heavy_check_mark:    |        [GPL 3.0](https://github.com/dmcable/spacexr/blob/master/LICENSE)        | Cable D.M. et al. Robust decomposition of cell type mixtures in spatial transcriptomics. Nature Biotechnology 40, 517–526 (2022). <https://doi.org/10.1038/s41587-021-00830-w>                                                      |
| [CARD](https://github.com/YingMa0107/CARD)                     |   :heavy_check_mark:    |      [GPL-3.0](https://github.com/YingMa0107/CARD/blob/master/LICENSE.md)       | Ma Y., Zhou X. Spatially informed cell type deconvolution for spatial transcriptomics. Nature Biotechnology 40, 1349–1359 (2022). <https://doi.org/10.1038/s41587-022-01273-7>                                                      |
| [DOT](https://github.com/saezlab/dot)                          |   :heavy_check_mark:    |         [GPL-3.0](https://github.com/saezlab/DOT/blob/main/LICENSE.md)          | Rahimi A. et al. DOT: a flexible multi-objective optimization framework for transferring features across single-cell and spatial omics. Nature Communications 15, 4994 (2024). <https://www.nature.com/articles/s41467-024-48868-z> |
