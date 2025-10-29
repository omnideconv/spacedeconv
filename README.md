# :rocket: spacedeconv <a href="https://omnideconv.github.io/spacedeconv"><img src="man/figures/logo.png" align="right" height="200" /></a>

[![R-CMD-check](https://github.com/omnideconv/spacedeconv/actions/workflows/test.yml/badge.svg)](https://github.com/omnideconv/spacedeconv/actions/workflows/test.yml)
[![docs](https://github.com/omnideconv/spacedeconv/actions/workflows/pkgdown.yml/badge.svg)](https://omnideconv.org/spacedeconv)
[![test-coverage](https://github.com/omnideconv/spacedeconv/actions/workflows/test-coverage.yml/badge.svg)](https://github.com/omnideconv/spacedeconv/actions/workflows/test-coverage.yml)
[![codecov](https://codecov.io/gh/omnideconv/spacedeconv/graph/badge.svg?token=OX9ZHSEP9L)](https://codecov.io/gh/omnideconv/spacedeconv)

*spacedeconv* is a unified interface to first- and second-generation deconvolution tools with focus on spatial transcriptomics datasets. The package is able to directly estimate celltype proportions of immune cells and can deconvolute any celltype if an annotated single-cell reference dataset is available.

## :arrow_down: Installation

**Note:** The current *spacedeconv* installation is only available for the **linux-64** platform.

Since many different packages need to be included, we highly recommend to install *spacedeconv* in a new Conda environment with the following commands.

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

Install the missing packages that are not available via conda as well as *spacedeconv*:
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

Single-cell data with cell-type annotation: _[SingleCellExperiment](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)_ (recommended), _[anndata](https://anndata.dynverse.org/)_ or _[Seurat](https://satijalab.org/seurat/)_

## :technologist: Usage

The main workflow consists of:

1. Reference signature computation using annotated single-cell data
2. Deconvolution
3. Visualization

### 1. Build a Signature Matrix

Build a cell type specific signature matrix from annotated single-cell reference data.

```r
signature <- spacedeconv::build_model(
  single_cell_object,
  cell_type_col = "celltype_major",
  method = "spotlight",
  assay_sc="cpm"
)
```

### 2. Deconvolution

While some methods are able to directly estimate immune cell abundances other tools require a custom reference signature computed in step 1).

```r
result <- spacedeconv::deconvolute(
  spatial_object,
  signature,
  method = "spotlight"
)
```

### 3. Visualization

*spacedeconv* includes highly-flexible visualization functions. A full explanation of all visualization options can be found in the visualization [vignette](articles/spacedeconv_visualization.html).

```r
plot_celltype(spe, cell_type="spotlight_B.cells")
```

## Available methods, Licenses, Citations

Note that, while *spacedeconv* itself is free ([GPL
3.0](https://github.com/omnideconv/omnideconv/blob/main/LICENSE)), you may
need to obtain a license to use the individual methods. See the table
below for more information. If you use this package in your work, please
cite both our package and the method(s) you are using.

> Constantin Zackl, Maria Zopoglou, Reto Stauffer, Markus Ausserhofer, Marieke E. Ijsselsteijn, Gregor Sturm, Noel Filipe da Cunha Carvalho de Miranda, Francesca Finotello. spacedeconv: deconvolution of tissue architecture from spatial transcriptomics, PREPRINT available at Research Square https://doi.org/10.21203/rs.3.rs-5102166/v1

| First-gen (immunedeconv)                                                                                                                                                               | First-gen mouse (immunedeconv)                                            | Second-gen (spatial methods)                                                                                                                                                                                                                                                           |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <ul><li>MPCcounter</li><li>EPIC</li><li>quanTIseq</li><li>xCell</li><li>CIBERSORT</li><li>CIBERSORT (abs.)</li><li>TIMER</li><li>ConsensusTME</li><li>ABIS</li><li>ESTIMATE</li> </ul> | <ul> <li>mMCPcounter</li><li>seqImmuCC</li><li>DCP</li><li>BASE</li></ul> | <ul><li>RCTD</li><li>SPOTlight</li><li>CARD</li><li>spatialDWLS</li><li>Cell2location</li><li>DOT</li></ul> |

# References

| Method                                                         |     signature      |                                     licence                                     | citation                                                                                                                                                                                                                                                          |
| -------------------------------------------------------------- | :----------------: | :-----------------------------------------------------------------------------: | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [immunedeconv](https://github.com/omnideconv/immunedeconv)     |        :x:         |      [BSD](https://github.com/omnideconv/immunedeconv/blob/master/LICENSE)      | Sturm, G., Finotello, F., Petitprez, F., Zhang, J. D., Baumbach, J., Fridman, W. H., ..., List, M., Aneichyk, T. (2019). Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. Bioinformatics, 35(14), i436-i445. |
| [spatialDWLS](https://github.com/RubD/Giotto/)                 | :heavy_check_mark: |            [MIT](https://github.com/RubD/Giotto/blob/master/LICENSE)            | Dong R, Yuan GC. SpatialDWLS: accurate deconvolution of spatial transcriptomic data. Genome Biol. 2021 May 10;22(1):145. doi: 10.1186/s13059-021-02362-7                                                                                                          |
| [cell2location](https://github.com/BayraktarLab/cell2location) | :heavy_check_mark: | [Apache-2.0](https://github.com/BayraktarLab/cell2location/blob/master/LICENSE) | Kleshchevnikov, V., Shmatko, A., Dann, E. et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nat Biotechnol (2022). <https://doi.org/10.1038/s41587-021-01139-4>                                                                       |
| [SPOTlight](https://github.com/MarcElosua/SPOTlight)           | :heavy_check_mark: |     [GPL 3.0](https://github.com/MarcElosua/SPOTlight/blob/main/LICENSE.md)     | Elosua-Bayes M, Nieto P, Mereu E, Gut I, Heyn H (2021): _SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes_. **Nucleic Acids Res** 49(9):e50. <doi:10.1093/nar/gkab043>.                              |
| [RCTD](https://github.com/dmcable/spacexr)                     | :heavy_check_mark: |        [GPL 3.0](https://github.com/dmcable/spacexr/blob/master/LICENSE)        | Cable, D.M., Murray, E., Zou, L.S. et al. Robust decomposition of cell type mixtures in spatial transcriptomics. Nat Biotechnol 40, 517–526 (2022). <https://doi.org/10.1038/s41587-021-00830-w>                                                                  |
| [CARD](https://github.com/YingMa0107/CARD)                     | :heavy_check_mark: |      [GPL-3.0](https://github.com/YingMa0107/CARD/blob/master/LICENSE.md)       | Ying Ma, and Xiang Zhou (2022). Spatially informed cell type deconvolution for spatial transcriptomics. <https://doi.org/10.1038/s41587-022-01273-7>                                                                                                              |
