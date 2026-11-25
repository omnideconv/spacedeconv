# :rocket: spacedeconv <a href="https://omnideconv.github.io/spacedeconv"><img src="man/figures/logo.png" align="right" height="200" /></a>

[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![R](https://img.shields.io/badge/R-4.3.3-276DC3?logo=r&logoColor=white)
![Python](https://img.shields.io/badge/Python-3.10-3776AB?logo=python&logoColor=white)
![Platform](https://img.shields.io/badge/platform-Linux%2064--bit-lightgrey)

_spacedeconv_ is a unified interface to first- and second-generation deconvolution tools with focus on spatial transcriptomics datasets. The package is able to directly estimate celltype proportions of immune cells and can deconvolute any celltype if an annotated single-cell reference dataset is available.

## :arrow_down: Installation

**Note:** The current _spacedeconv_ installation is only available for the **linux-64** platform.

Since many different packages need to be included, we highly recommend to install _spacedeconv_ in a new Conda environment with the following commands.

First, a tool for fast dependency resolution is needed, therefore we recommend installing [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) (or mamba / conda with libmamba solver) if not already available:

```r
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

Download the environment.yml file of this github repo:

```r
wget https://raw.githubusercontent.com/omnideconv/spacedeconv/main/environment.yml -O environment.yml
```

Create a new environment called "spacedeconv-env" via micromamba with the environment.yml file and flexible channel priority:

```r
micromamba env create -f environment.yml --channel-priority flexible
```

Start R inside the spacedeconv-env conda environment:

```r
conda activate spacedeconv-env
R
```

Install the missing packages that are not available via conda as well as _spacedeconv_:

```r
pak::pkg_install("cran::wrMisc", dependencies = FALSE, upgrade = FALSE)
pak::pkg_install("drieslab/Giotto@v3.3.2", dependencies = FALSE, upgrade = FALSE)
pak::pkg_install("xuranw/MuSiC@f21fe67f5670d5e9fca0ad7550abaae3423eb59c", dependencies = FALSE, upgrade = FALSE)
pak::pkg_install("YingMa0107/CARD@2d64b91abb5cdd0c7f576b1c5d4727c84e7c93a0", dependencies = FALSE, upgrade = FALSE)
pak::pkg_install("saezlab/DOT", dependencies = FALSE, upgrade = FALSE)
pak::pkg_install("federicomarini/quantiseqr@1bc02b04350d190ed29de3f9c179f92ad70337af", dependencies = FALSE, upgrade = FALSE)
pak::pkg_install("saezlab/OmnipathR@v4.0.4", dependencies = FALSE, upgrade = FALSE)
pak::pkg_install("omnideconv/spacedeconv", dependencies = FALSE, upgrade = FALSE)
```

## :sparkles: Features

- Unified access to first- and second-generation deconvolution tools
- Direct deconvolution of immune cells
- Compute custom reference signatures to deconvolute any celltype
- Flexible visualization functions
- Resource optimization
- Pathway and Transcription Factor analysis (decoupleR integration)
- Ligand-Receptor analysis
- Easy integration into spatial transcriptomics workflows

<div align="center">
  <img src="man/figures/spacedeconv_summary.png" alt="spacedeconv summary diagram" width="800"/>
</div>

<div align="center">
  <img src="man/figures/spacedeconv_summary.png" alt="spacedeconv summary diagram" width="800"/>
</div>

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

Alternatively, load data from a repository such as Zenodo. Spatial data should be provided as a `SpatialExperiment` and single-cell reference data as a `SingleCellExperiment`. The example below downloads the same `spatial_data_3` and `single_cell_data_3` objects from [Zenodo](https://doi.org/10.5281/zenodo.24108820) instead of loading them with `data()`.

```r
download.file(
  "https://zenodo.org/records/24108820/files/spatial_data_3.RData?download=1",
  destfile = "spatial_data_3.RData",
  mode = "wb"
)
load("spatial_data_3.RData")

download.file(
  "https://zenodo.org/records/24108820/files/single_cell_data_3.RData?download=1",
  destfile = "single_cell_data_3.RData",
  mode = "wb"
)
load("single_cell_data_3.RData")
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

You can see all supported deconvolution methods with `spacedeconv::deconvolution_methods`.

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
  density = FALSE
)
```

## Available methods, licenses, and citations

We focused on supporting methods that have demonstrated robust performance
in published benchmarking analyses, while also considering the feasibility
of integration, installation, and maintenance. Results from external
deconvolution tools not directly supported by _spacedeconv_ can also be imported
for downstream analysis and visualization, as demonstrated in the
[external deconvolution results vignette](https://felixpetschko.github.io/spacedeconv/articles/spacedeconv_import_deconvolution_results.html).

To balance a minimal installation burden for users with convenient access
to high-performing methods, we plan to consider incorporating newly published
methods upon request or through contributions from their developers,
provided that independent or published benchmarking evidence demonstrates a substantial performance gain over the methods already supported.

Note that, while _spacedeconv_ itself is free ([GPL
3.0](https://github.com/omnideconv/spacedeconv/blob/main/LICENSE)), you may
need to obtain a license to use the individual methods. See the table
below for more information. If you use this package in your work, please
cite both our package and the method(s) you are using.

> Constantin Zackl, Maria Zopoglou, Reto Stauffer, Markus Ausserhofer, Marieke E. Ijsselsteijn, Gregor Sturm, Noel Filipe da Cunha Carvalho de Miranda, Francesca Finotello. spacedeconv: deconvolution of tissue architecture from spatial transcriptomics, PREPRINT available at Research Square https://doi.org/10.21203/rs.3.rs-5102166/v1.

First-generation methods use predefined signatures or gene sets and are accessed through immunedeconv ([Sturm et al., 2019](https://doi.org/10.1093/bioinformatics/btz363); [Merotto et al., 2024](https://doi.org/10.1093/bioadv/vbae032) for mouse deconvolution). Please cite both immunedeconv and the original method. Second-generation methods use a user-provided annotated single-cell or single-nucleus reference.

“Originally developed for” indicates the method’s original application. Bulk methods are applied to spatial spots as individual samples in spacedeconv. Species entries describe the standard references available through spacedeconv. For reference-dependent methods, species support depends on a suitable annotated reference and compatible gene identifiers.

All methods can run on CPUs; for cell2location, GPU acceleration is supported and recommended. Memory requirements depend on the method and dataset size, including the numbers of spots, genes, and reference cells.

| Method | Generation | Originally developed for | Species | License | Citation |
|---|---|---|---|---|---|
| [MCP-counter](https://github.com/ebecht/MCPcounter) | First | Bulk | Human | [GPL-3.0](https://github.com/ebecht/MCPcounter/blob/master/Source/License) | Becht E. et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology 17, 218 (2016). <https://doi.org/10.1186/s13059-016-1070-5> |
| [EPIC](https://github.com/GfellerLab/EPIC) | First | Bulk | Human | [Academic non-commercial research](https://github.com/GfellerLab/EPIC/blob/master/LICENSE) | Racle J. et al. Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. eLife 6, e26476 (2017). <https://doi.org/10.7554/eLife.26476> |
| [quanTIseq](https://icbi.i-med.ac.at/software/quantiseq/doc/index.html) | First | Bulk | Human | [BSD-3-Clause](https://icbi.i-med.ac.at/software/quantiseq/doc/index.html) | Finotello F. et al. Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome Medicine 11, 34 (2019). <https://doi.org/10.1186/s13073-019-0638-6> |
| [xCell](https://github.com/dviraran/xCell) | First | Bulk | Human | [GPL-3.0](https://github.com/dviraran/xCell/blob/master/DESCRIPTION) | Aran D. et al. xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biology 18, 220 (2017). <https://doi.org/10.1186/s13059-017-1349-1> |
| [CIBERSORT](https://cibersortx.stanford.edu/) | First | Bulk | Human (LM22) | [Non-commercial terms](https://cibersortx.stanford.edu/) | Newman A.M. et al. Robust enumeration of cell subsets from tissue expression profiles. Nature Methods 12, 453–457 (2015). <https://doi.org/10.1038/nmeth.3337> |
| [TIMER](https://doi.org/10.1186/s13059-016-1028-7) | First | Bulk | Human | [GPL-2.0](https://omnideconv.org/immunedeconv/#available-methods-licenses-citations) | Li B. et al. Comprehensive analyses of tumor immunity: implications for cancer immunotherapy. Genome Biology 17, 174 (2016). <https://doi.org/10.1186/s13059-016-1028-7> |
| [ConsensusTME](https://github.com/cansysbio/ConsensusTME) | First | Bulk | Human | [GPL-3.0](https://github.com/cansysbio/ConsensusTME/blob/master/LICENSE.md) | Jiménez-Sánchez A. et al. Comprehensive Benchmarking and Integration of Tumor Microenvironment Cell Estimation Methods. Cancer Research 79(24), 6238–6246 (2019). <https://doi.org/10.1158/0008-5472.CAN-18-3560> |
| [ABIS](https://github.com/giannimonaco/ABIS) | First | Bulk | Human | [GPL-2.0](https://github.com/giannimonaco/ABIS) | Monaco G. et al. RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types. Cell Reports 26(6), 1627–1640.e7 (2019). <https://doi.org/10.1016/j.celrep.2019.01.041> |
| [ESTIMATE](https://bioinformatics.mdanderson.org/public-software/estimate/) | First | Bulk | Human | [GPL-2.0](https://bioinformatics.mdanderson.org/public-software/estimate/) | Yoshihara K. et al. Inferring tumour purity and stromal and immune cell admixture from expression data. Nature Communications 4, 2612 (2013). <https://doi.org/10.1038/ncomms3612> |
| [mMCP-counter](https://github.com/cit-bioinfo/mMCP-counter) | First | Bulk | Mouse | [GPL-3.0](https://github.com/cit-bioinfo/mMCP-counter/blob/master/LICENSE.md) | Petitprez F. et al. The murine Microenvironment Cell Population counter method to estimate abundance of tissue-infiltrating immune and stromal cell populations in murine samples using gene expression. Genome Medicine 12, 86 (2020). <https://doi.org/10.1186/s13073-020-00783-w> |
| [seqImmuCC](https://doi.org/10.3389/fimmu.2018.01286) | First | Bulk | Mouse | [Non-commercial terms](https://omnideconv.org/immunedeconv/#available-methods-licenses-citations) | Chen Z. et al. seq-ImmuCC: Cell-Centric View of Tissue Transcriptome Measuring Cellular Compositions of Immune Microenvironment From Mouse RNA-Seq Data. Frontiers in Immunology 9, 1286 (2018). <https://doi.org/10.3389/fimmu.2018.01286> |
| [DCQ](https://cran.r-project.org/package=ComICS) | First | Bulk | Mouse | [GPL-2.0](https://cran.r-project.org/package=ComICS) | Altboum Z. et al. Digital cell quantification identifies global immune cell dynamics during influenza infection. Molecular Systems Biology 10, 720 (2014). <https://doi.org/10.1002/msb.134947> |
| [BASE](https://doi.org/10.1038/ncomms10248) | First | Bulk | Mouse | [See upstream notes](https://omnideconv.org/immunedeconv/#available-methods-licenses-citations) | Varn F.S. et al. Integrative analysis of breast cancer reveals prognostic haematopoietic activity and patient-specific immune response profiles. Nature Communications 7, 10248 (2016). <https://doi.org/10.1038/ncomms10248> |
| [spatialDWLS](https://github.com/giotto-suite/Giotto/tree/v3.3.2) | Second | Spatial | Reference-dependent | [GPL-3.0-or-later](https://github.com/giotto-suite/Giotto/blob/v3.3.2/DESCRIPTION) | Dong R., Yuan G.C. SpatialDWLS: accurate deconvolution of spatial transcriptomic data. Genome Biology 22, 145 (2021). <https://doi.org/10.1186/s13059-021-02362-7> |
| [cell2location](https://github.com/BayraktarLab/cell2location) | Second | Spatial | Reference-dependent | [Apache-2.0](https://github.com/BayraktarLab/cell2location/blob/master/LICENSE) | Kleshchevnikov V. et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nature Biotechnology 40, 661–671 (2022). <https://doi.org/10.1038/s41587-021-01139-4> |
| [SPOTlight](https://github.com/MarcElosua/SPOTlight) | Second | Spatial | Reference-dependent | [GPL 3.0](https://github.com/MarcElosua/SPOTlight/blob/devel/LICENSE.md) | Elosua-Bayes M. et al. SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes. Nucleic Acids Research 49(9):e50 (2021). <https://doi.org/10.1093/nar/gkab043>. |
| [RCTD](https://github.com/dmcable/spacexr) | Second | Spatial | Reference-dependent | [GPL 3.0](https://github.com/dmcable/spacexr/blob/master/LICENSE) | Cable D.M. et al. Robust decomposition of cell type mixtures in spatial transcriptomics. Nature Biotechnology 40, 517–526 (2022). <https://doi.org/10.1038/s41587-021-00830-w> |
| [CARD](https://github.com/YMa-lab/CARD) | Second | Spatial | Reference-dependent | [GPL-3.0-or-later](https://github.com/YMa-lab/CARD/blob/master/DESCRIPTION) | Ma Y., Zhou X. Spatially informed cell type deconvolution for spatial transcriptomics. Nature Biotechnology 40, 1349–1359 (2022). <https://doi.org/10.1038/s41587-022-01273-7> |
| [DOT](https://github.com/saezlab/dot) | Second | Spatial | Reference-dependent | [GPL-3.0-or-later](https://github.com/saezlab/DOT/blob/main/DESCRIPTION) | Rahimi A. et al. DOT: a flexible multi-objective optimization framework for transferring features across single-cell and spatial omics. Nature Communications 15, 4994 (2024). <https://www.nature.com/articles/s41467-024-48868-z> |
| [Rectangle](https://github.com/ComputationalBiomedicineGroup/Rectangle) | Second | Bulk / spatial | Reference-dependent | [BSD-3-Clause OR Commercial](https://github.com/ComputationalBiomedicineGroup/Rectangle/blob/main/LICENSE) | Eder B, Rigato I, Dietrich A, Merotto L, Sturm G, Treis T, List M, Theis FJ, Finotello F. Rectangle: robust and scalable multiscale deconvolution informed by single-cell RNA sequencing data. bioRxiv. 2026. <https://doi.org/10.64898/2026.07.07.736950> |
