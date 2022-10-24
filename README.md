# SpaceDeconv

[![R-CMD-check](https://github.com/omnideconv/SpaceDeconv/actions/workflows/test.yml/badge.svg)](https://github.com/omnideconv/SpaceDeconv/actions/workflows/test.yml)
[![docs](https://github.com/omnideconv/SpaceDeconv/actions/workflows/pkgdown.yml/badge.svg)](https://github.com/omnideconv/SpaceDeconv/actions/workflows/pkgdown.yml)

SpaceDeconv is a unified interface to 31 deconvolution tools with focus on spatial transcriptomics datasets. In total 17 second-generation deconvolution tools are included, enabling deconvolution of any cell types when single-cell reference data is available. Additionally 10 first-generation tools, which are focusing on deconvolution of immune cells, are available as well as 4 first-generation methods optimised for mouse data.

## Installation

There are two ways to install `SpaceDeconv`:

- The _minimal_ installation installs only the dependencies required for the basic functionalities. All deconvolution methods need to be installed on-demand.
- The _complete_ installation installs all dependencies including all deconvolution methods. This may take a considerable time.

Since not all dependencies are on CRAN or Bioconductor, `SpaceDeconv`is available from GitHub only. We recommend installing trough the pak package manager:

```r
# install the pak package manager
install.packages("pak")

# minimal installation
pak::pkg_install("omnideconv/SpaceDeconv")

# complete installation, including Python dependencies
pak::pkg_install("omnideconv/SpaceDeconv", dependencies=TRUE)
# SpaceDeconv::install_all_python()
```

## Usage

SpaceDeconv offers convenient access to perform first- and second-generation deconvolution on spatial transcriptomics datasets. While deconvolution can be performed directly with first-generation methods, second-generation algorithms require an additional annotated single-cell reference. A full list of deconvolution tools can be accessed by `SpaceDeconv::deconvolution_algorithms` or in the [FAQ](articles/SpaceDeconv_faq.html).

### Data requirements

- _[SpatialExperiment](https://bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html)_, will be deconvoluted
- _[SingleCellExperiment](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)_ (recommended), _[anndata](https://anndata.dynverse.org/)_ or _[Seurat](https://satijalab.org/seurat/)_ containing cell type information

The main functions of SpaceDeconv are used to build a signature matrix from annotated single-cell transcriptomics data and deconvolute spatially resolved transcriptomics datasets. The basic workflow consists of:

### 1. Build a Signature Matrix

Build a cell type specific signature matrix from annotated single-cell reference data.

```r
signature <- SpaceDeconv::build_model(
  single_cell_object,
  cell_type_col = "celltype_major",
  method = "spotlight",
  assay_sc="cpm"
)
```

### 2. Deconvolution

To perform a deconvolution a SpatialExperiment object is required. Some methods additionally require a cell-type specific reference signature which can be calculated by `SpaceDeconv::build_model()`. By default the deconvolution results are added to the SpatialExperiment object to simplify the visualization. You can obtain the results in table form by setting `return_object=FALSE`.

```r
# save the results to an annotated SpatialExperiment
result <- SpaceDeconv::deconvolute(
  spatial_object,
  signature,
  method = "spotlight"
)
```

### 3. Visualization

SpaceDeconv includes multiple visualization functions.

```r
# sample does refer to the first column of ColData(spe)
# for cell_type input a celltype present in the deconvolution result
plot_celltype(spe, cell_type="spotlight_B.cells")
```

## Additional Requirements

Most methods do not require additional software/tokens, but there are a few exceptions:

- A working version of Docker is required for CIBERSORTx
- A token for CIBERSORTx is required from this website:
  <https://cibersortx.stanford.edu/>
- The CIBERSORT source code is required for BSeq-sc (see tutorial in
  ?omnideconv::bseqsc_config)
- SpatialExperiment requires `magick` to be installed: `sudo apt-get install libmagick++-dev`

## Available methods, Licenses, Citations

Note that, while _SpaceDeconv_ itself is free ([GPL
3.0](https://github.com/omnideconv/omnideconv/blob/main/LICENSE)), you may
need to obtain a license to use the individual methods. See the table
below for more information. If you use this package in your work, please
cite both our package and the method(s) you are using.

| First-gen (immunedeconv)                                                                                                                                                               | First-gen mouse (immunedeconv)                                            | Second-gen (omnideconv + spatial Methods)                                                                                                                                                                                                                                                           |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <ul><li>MPCcounter</li><li>EPIC</li><li>quanTIseq</li><li>xCell</li><li>CIBERSORT</li><li>CIBERSORT (abs.)</li><li>TIMER</li><li>ConsensusTME</li><li>ABIS</li><li>ESTIMATE</li> </ul> | <ul> <li>mMCPcounter</li><li>seqImmuCC</li><li>DCP</li><li>BASE</li></ul> | <ul><li>RCTD</li><li>SPOTlight</li><li>CARD</li><li>spatialDWLS</li><li>Cell2location</li><li>AutoGeneS</li><li>BayesPrism</li><li>Bisque</li><li>Bisque</li><li>Bseq-sc</li><li>CIBERSORTx</li><li>CDseq</li><li>CPM</li><li>DWLS</li><li>MOMF</li><li>MuSiC</li><li>Scaden</li><li>SCDC</li></ul> |

# References

| Method                                                         |     signature      |     organism     |                                     licence                                     | citation                                                                                                                                                                                                                                                          |
| -------------------------------------------------------------- | :----------------: | :--------------: | :-----------------------------------------------------------------------------: | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [omnideconv](https://github.com/omnideconv/omnideconv)         | :heavy_check_mark: |      human       |       [GPL-3](https://github.com/omnideconv/omnideconv/blob/main/LICENSE)       | Citation will follow                                                                                                                                                                                                                                              |
| [immunedeconv](https://github.com/omnideconv/immunedeconv)     |        :x:         | humand and mouse |      [BSD](https://github.com/omnideconv/immunedeconv/blob/master/LICENSE)      | Sturm, G., Finotello, F., Petitprez, F., Zhang, J. D., Baumbach, J., Fridman, W. H., ..., List, M., Aneichyk, T. (2019). Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. Bioinformatics, 35(14), i436-i445. |
| [spatialDWLS](https://github.com/RubD/Giotto/)                 | :heavy_check_mark: |      human       |            [MIT](https://github.com/RubD/Giotto/blob/master/LICENSE)            | Dong R, Yuan GC. SpatialDWLS: accurate deconvolution of spatial transcriptomic data. Genome Biol. 2021 May 10;22(1):145. doi: 10.1186/s13059-021-02362-7                                                                                                          |
| [cell2location](https://github.com/BayraktarLab/cell2location) | :heavy_check_mark: |      human       | [Apache-2.0](https://github.com/BayraktarLab/cell2location/blob/master/LICENSE) | Kleshchevnikov, V., Shmatko, A., Dann, E. et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nat Biotechnol (2022). <https://doi.org/10.1038/s41587-021-01139-4>                                                                       |
| [SPOTlight](https://github.com/MarcElosua/SPOTlight)           | :heavy_check_mark: |      human       |     [GPL 3.0](https://github.com/MarcElosua/SPOTlight/blob/main/LICENSE.md)     | Elosua-Bayes M, Nieto P, Mereu E, Gut I, Heyn H (2021): _SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes_. **Nucleic Acids Res** 49(9):e50. <doi:10.1093/nar/gkab043>.                              |
| [RCTD](https://github.com/dmcable/spacexr)                     | :heavy_check_mark: |      human       |        [GPL 3.0](https://github.com/dmcable/spacexr/blob/master/LICENSE)        | Cable, D.M., Murray, E., Zou, L.S. et al. Robust decomposition of cell type mixtures in spatial transcriptomics. Nat Biotechnol 40, 517–526 (2022). <https://doi.org/10.1038/s41587-021-00830-w>                                                                  |
| [CARD](https://github.com/YingMa0107/CARD)                     | :heavy_check_mark: |      human       |      [GPL-3.0](https://github.com/YingMa0107/CARD/blob/master/LICENSE.md)       | Ying Ma, and Xiang Zhou (2022). Spatially informed cell type deconvolution for spatial transcriptomics. <https://doi.org/10.1038/s41587-022-01273-7>                                                                                                              |
