# SpaceDeconv

SpaceDeconv is a unified interface to 31 deconvolution tools with focus on spatial transcriptomics datasets. In total 17 second-generation deconvolution tools are available, enabling deconvolution of any cell types when single-cell reference data is available. Additionally 10 first-generation tools, which are focussing on deconvolution of immune cells, are available as well as 4 first-generation methods optimised for mouse data.

## Installation 

There are two ways to install `SpaceDeconv`:
* The _minimal_ installation installs only the dependencies required for the basic functionalities. All deconvolution methods need to be installed on-demand. 
* The _complete_ installtation installs all dependencies including all deconvolution methods. This may take a considerable time. 

Since not all dependencies are on CRAN or Bioconductor, `SpaceDeconv`is available from GitHub only. We recommend installing trough the pak package manager: 

``` r
# install the pak package manager
install.packages("pak")

# minimal installation 
pak::pkg_install("omnideconv/SpaceDeconv")

# complete installation, including Python dependencies
pak::pkg_install("omnideconv/SpaceDeconv", dependencies=TRUE)
# SpaceDeconv::install_all_python()
``` 

## Usage 

### Data requirements
SpaceDeconv offers convenient access to perform first- and second-generation deconvolution on spatial transcriptomics datasets. 

* _SpatialExperiment_, will be deconvoluted
* _SingleCellExperiment_ containing cell type information, for second-generation tools

### Load Spatial Dataset
You can load 10X Visium Data by providing the spaceranger output folder. It is further possible to run SpaceDeconv with manually created SpatialExperiments. See the SpatialExperiment [Documentation](https://github.com/drighelli/SpatialExperiment) for further details. 
``` r
spe <- SpatialExperiment::read10xVisium("path_to_spaceranger_output")
``` 

### Normalization 
SpaceDeconv offers an additional function for convenient normalization of SpatialExperiments. The normalization is saved in a new assay, so make sure the correct data is used during deconvolution by providing the desired assay with the parameters `assay_sc` and `assay_sp`.
``` r
spe <- SpaceDeconv::normalize(spe, method="cpm")

```

### Build a Signature Matrix 
``` r
signature <- SpaceDeconv::build_model(
  single_cell_object, 
  cell_type_col = "cell_ontology_class", 
  method = "spotlight", 
  assay_sc="cpm"
)
```
### Deconvolution
To perform a deconvolution a SpatialExperiment object is required. Some methods additionally require a cell-type specific reference signature which can be calculated by `SpaceDeconv::build_model()`. By default the deconvolution results are added to the SpatialExperiment object to simplify the visualization. You can obtain the results in table form by setting `return_object=FALSE`.
```r
# save the results to an annotated SpatialExperiment
result <- SpaceDeconv::deconvolute(
  spatial_object, 
  signature, 
  method = "spotlight"
)

# return deconvolution results as a table
result <- SpaceDeconv::deconvolute(
  spatial_object, 
  signature, 
  method = "spotlight", 
  return_object = FALSE
)
```


### Visualization 
SpaceDeconv includes multiple visualization functions. 
```r
# sample does refer to the first column of ColData(spe)
# for cell_type input a celltype present in the deconvolution result
plot_celltype(spe, sample="sample01", cell_type="B.cells")

# threshold changes the minimum cell type fraction 
# for a cell to be considered present in a specific spot
plot_cells_per_spot(spe, plot_type = "spatial", threshold=0.01)
```

## Additional Requirements
Most methods do not require additional software/tokens, but there are a feq exceptions: 

- A working version of Docker is required for CIBERSORTx
- A token for CIBERSORTx is required from this website:
  <https://cibersortx.stanford.edu/>
- The CIBERSORT source code is required for BSeq-sc (see tutorial in
  ?omnideconv::bseqsc_config)
  
## Available methods, Licenses, Citations

Note that, while _SpaceDeconv_ itself is free ([GPL
3.0](https://github.com/omnideconv/omnideconv/blob/main/LICENSE)), you may
need to obtain a license to use the individual methods. See the table
below for more information. If you use this package in your work, please
cite both our package and the method(s) you are using.


| Method | reference | organism | licence | citation |
| --- | :---: | :---: | :---: | --- |
| [quanTIseq](http://icbi.at/software/quantiseq/doc/index.html)| :x: | human | free ([BSD](https://github.com/omnideconv/immunedeconv/blob/master/LICENSE.md)) | Finotello, F., Mayer, C., Plattner, C., Laschober, G., Rieder, D., Hackl, H., ..., Sopper, S. (2019). Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome medicine, 11(1), 34. https://doi.org/10.1186/s13073-019-0638-6 |
| [TIMER](http://cistrome.org/TIMER/)| :x: | human | free ([GPL 2.0](http://cistrome.org/TIMER/download.html)) | Li, B., Severson, E., Pignon, J.-C., Zhao, H., Li, T., Novak, J., … Liu, X. S. (2016). Comprehensive analyses of tumor immunity: implications for cancer immunotherapy. Genome Biology, 17(1), 174. https://doi.org/10.1186/s13059-016-1028-7 |
| [CIBERSORT](https://cibersort.stanford.edu/) | :x: | human    | free for non-commerical use only | Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., … Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453–457. https://doi.org/10.1038/nmeth.3337 |
| [MCPCounter](https://github.com/ebecht/MCPcounter) | :x: | human    | free ([GPL 3.0](https://github.com/ebecht/MCPcounter/blob/master/Source/License)) | Becht, E., Giraldo, N. A., Lacroix, L., Buttard, B., Elarouci, N., Petitprez, F., … de Reyniès, A. (2016). Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology, 17(1), 218. https://doi.org/10.1186/s13059-016-1070-5 |
| [xCell](http://xcell.ucsf.edu/) | :x: | human    | free ([GPL 3.0](https://github.com/dviraran/xCell/blob/master/DESCRIPTION))   | Aran, D., Hu, Z., & Butte, A. J. (2017). xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biology, 18(1), 220. https://doi.org/10.1186/s13059-017-1349-1|
| [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/) | :x: | human    | free for non-commercial use only ([Academic License](https://github.com/GfellerLab/EPIC/blob/master/LICENSE)) | Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E., & Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. ELife, 6, e26476. https://doi.org/10.7554/eLife.26476 |
| [ESTIMATE](https://gfellerlab.shinyapps.io/EPIC_1-1/) | :x: | human  | free ([GPL 2.0](https://bioinformatics.mdanderson.org/public-software/estimate/))   | Yoshihara, K., Shahmoradgoli, M., Martínez, E., Vegesna, R., Kim, H., Torres-Garcia, W., Treviño, V., Shen, H., Laird, P. W., Levine, D. A., Carter, S. L., Getz, G., Stemke-Hale, K., Mills, G. B., & Verhaak, R. G. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. Nature communications, 4, 2612. https://doi.org/10.1038/ncomms3612 |
| [ABIS](https://giannimonaco.shinyapps.io/ABIS/) | :x: | human    | free ([GPL 2.0](https://github.com/giannimonaco/ABIS)) | Monaco, G., Lee, B., Xu, W., Mustafah, S., Hwang, Y. Y., ..., Larbi, A. (2019). RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types. Cell reports, 26(6), 1627–1640.e7. https://doi.org/10.1016/j.celrep.2019.01.041                                                                                                             |
| [ConsensusTME](https://olliecast.shinyapps.io/Deconvolution_Benchmarking/) | :x: |  human    | free ([GPL 3.0](https://github.com/cansysbio/ConsensusTME/blob/master/LICENSE.md))                            | Jiménez-Sánchez, A., Cast, O., & Miller, M. L. (2019). Comprehensive Benchmarking and Integration of Tumor Microenvironment Cell Estimation Methods. Cancer research, 79(24), 6238–6246. https://doi.org/10.1158/0008-5472.CAN-18-3560  |
| [mMCPCounter](https://github.com/cit-bioinfo/mMCP-counter) | :x:  | mouse    | free ([GPL 3.0](https://github.com/cit-bioinfo/mMCP-counter/blob/master/LICENSE.md))  | Petitprez, F., Levy, S., Sun, C. M., Meylan, M., ..., de Reyniès, A. (2020). The murine Microenvironment Cell Population counter method to estimate abundance of tissue-infiltrating immune and stromal cell populations in murine samples using gene expression. Genome medicine, 12(1), 86. https://doi.org/10.1186/s13073-020-00783-w  |
| [seqImmuCC](218.4.234.74:3200/immune/) | :x: | mouse    | free for non-commerical use only  | Chen, Z., Quan, L., Huang, A., Zhao, Q., Yuan, Y., Yuan, X., ..., Wu, A. (2018). seq-ImmuCC: Cell-Centric View of Tissue Transcriptome Measuring Cellular Compositions of Immune Microenvironment From Mouse RNA-Seq Data. Frontiers in immunology, 9, 1286. https://doi.org/10.3389/fimmu.2018.01286                                                                                |
| [DCQ](http://dcq.tau.ac.il/)  | :x: | mouse    | free ([GPL 2.0](https://cran.r-project.org/web/packages/ComICS/index.html))                                   | Altboum, Z., Steuerman, Y., David, E., Barnett-Itzhaki, Z., Valadarsky, L., ..., Amit, I. (2014). Digital cell quantification identifies global immune cell dynamics during influenza infection. Molecular systems biology, 10(2), 720. https://doi.org/10.1002/msb.134947  |
| BASE  | :x:  | mouse    | free                                                                                                          | Varn, F. S., Andrews, E. H., Mullins, D. W., & Cheng, C. (2016). Integrative analysis of breast cancer reveals prognostic haematopoietic activity and patient-specific immune response profiles. Nature communications, 7, 10248. https://doi.org/10.1038/ncomms10248  |
| [AutoGeneS](https://github.com/theislab/AutoGeneS/) | :heavy_check_mark: | human | free ([MIT](https://github.com/theislab/AutoGeneS/blob/master/LICENSE))             | Aliee, H., & Theis, F. (2021). AutoGeneS: Automatic gene selection using multi-objective optimization for RNA-seq deconvolution. <https://doi.org/10.1101/2020.02.21.940650> |
| [Bisque](https://github.com/cozygene/bisque)       | :heavy_check_mark: | human    | free ([GPL 3.0](https://github.com/cozygene/bisque/blob/master/DESCRIPTION))        | Jew, B., Alvarez, M., Rahmani, E., Miao, Z., Ko, A., Garske, K. M., Sul, J. H., Pietiläinen, K. H., Pajukanta, P., & Halperin, E. (2020). Publisher Correction: Accurate estimation of cell composition in bulk expression through robust integration of single-cell information. Nature Communications, 11(1), 2891. <https://doi.org/10.1038/s41467-020-16607-9>                            |
| [BSeq-sc](https://github.com/shenorrLab/bseqsc)    | :heavy_check_mark: | human    | free ([GPL 2.0](https://github.com/shenorrLab/bseqsc/blob/master/DESCRIPTION))      | Baron, M., Veres, A., Wolock, S. L., Faust, A. L., Gaujoux, R., Vetere, A., Ryu, J. H., Wagner, B. K., Shen-Orr, S. S., Klein, A. M., Melton, D. A., & Yanai, I. (2016). A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. In Cell Systems (Vol. 3, Issue 4, pp. 346–360.e4). <https://doi.org/10.1016/j.cels.2016.08.011> |
| [CDSeq](https://github.com/kkang7/CDSeq_R_Package)  | :heavy_check_mark: | human   | free ([GPL 3.0](https://github.com/kkang7/CDSeq_R_Package/blob/master/DESCRIPTION)) | Kang, K., Huang, C., Li, Y. et al. CDSeqR: fast complete deconvolution for gene expression data from bulk tissues. BMC Bioinformatics 22, 262 (2021). <https://doi.org/10.1186/s12859-021-04186-5>                                                                                                                                                                                            |
| [CIBERSORTx](https://cibersortx.stanford.edu/)   | :heavy_check_mark: | human      | free for non-commerical use only                                                    | Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., Hoang, C. D., Diehn, M., & Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453–457. <https://doi.org/10.1038/nmeth.3337>                                                                                                                        |
| [CPM](https://github.com/amitfrish/scBio)    | :heavy_check_mark: | human          | free ([GPL 2.0](https://github.com/amitfrish/scBio/blob/master/DESCRIPTION))        | Frishberg, A., Peshes-Yaloz, N., Cohn, O., Rosentul, D., Steuerman, Y., Valadarsky, L., Yankovitz, G., Mandelboim, M., Iraqi, F. A., Amit, I., Mayo, L., Bacharach, E., & Gat-Viks, I. (2019). Cell composition analysis of bulk genomics using single-cell data. Nature Methods, 16(4), 327–332. <https://doi.org/10.1038/s41592-019-0355-5>    | :heavy_check_mark: | human                                             |
| [DWLS](https://bitbucket.org/yuanlab/dwls/src/master/) | :heavy_check_mark: | human | free ([GPL](https://bitbucket.org/yuanlab/dwls/src/master/DESCRIPTION))             | Tsoucas, D., Dong, R., Chen, H., Zhu, Q., Guo, G., & Yuan, G.-C. (2019). Accurate estimation of cell-type composition from gene expression data. Nature Communications, 10(1), 2975. <https://doi.org/10.1038/s41467-019-10802-z>                                                                                                                                                             |
| [MOMF](https://github.com/sqsun/MOMF)     | :heavy_check_mark: | human             | free ([GPL 3.0](https://github.com/sqsun/MOMF/blob/master/LICENSE.md))              | Xifang Sun, Shiquan Sun, and Sheng Yang. An efficient and flexible method for deconvoluting bulk RNAseq data with single-cell RNAseq data, 2019, DIO: 10.5281/zenodo.3373980                                                                                                                                                                                                                  |
| [MuSiC](https://github.com/xuranw/MuSiC/)   | :heavy_check_mark: | human           | free ([GPL 3.0](https://github.com/xuranw/MuSiC/blob/master/LICENSE))               | Wang, X., Park, J., Susztak, K., Zhang, N. R., & Li, M. (2019). Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nature Communications, 10(1), 380. <https://doi.org/10.1038/s41467-018-08023-x>                                                                                                                                                      |
| [Scaden](https://github.com/KevinMenden/scaden)    | :heavy_check_mark: | human    | free ([MIT](https://github.com/KevinMenden/scaden/blob/master/LICENSE))             | Menden, K., Marouf, M., Oller, S., Dalmia, A., Kloiber, K., Heutink, P., & Bonn, S. (n.d.). Deep-learning-based cell composition analysis from tissue expression profiles. <https://doi.org/10.1101/659227>                                                                                                                                                                                   |
| [SCDC](https://github.com/meichendong/SCDC)   | :heavy_check_mark: | human         | ([MIT](https://github.com/meichendong/SCDC/blob/master/README.md))                  | Dong, M., Thennavan, A., Urrutia, E., Li, Y., Perou, C. M., Zou, F., & Jiang, Y. (2020). SCDC: bulk gene expression deconvolution by multiple single-cell RNA sequencing references. Briefings in Bioinformatics. <https://doi.org/10.1093/bib/bbz166> |
| [spatialDWLS](https://github.com/RubD/Giotto/) | :heavy_check_mark: | human | [MIT](https://github.com/RubD/Giotto/blob/master/LICENSE) | Dong R, Yuan GC. SpatialDWLS: accurate deconvolution of spatial transcriptomic data. Genome Biol. 2021 May 10;22(1):145. doi: 10.1186/s13059-021-02362-7 |
| [cell2location](https://github.com/BayraktarLab/cell2location) | :heavy_check_mark: | human | [Apache-2.0](https://github.com/BayraktarLab/cell2location/blob/master/LICENSE) | Kleshchevnikov, V., Shmatko, A., Dann, E. et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nat Biotechnol (2022). <https://doi.org/10.1038/s41587-021-01139-4>|
| [SPOTlight](https://github.com/MarcElosua/SPOTlight) | :heavy_check_mark: | human | [GPL 3.0](https://github.com/MarcElosua/SPOTlight/blob/main/LICENSE.md) | Elosua-Bayes M, Nieto P, Mereu E, Gut I, Heyn H (2021): *SPOTlight: seeded NMF regression to deconvolute spatial transcriptomics spots with single-cell transcriptomes*. **Nucleic Acids Res** 49(9):e50. <doi:10.1093/nar/gkab043>. |
| [RCTD](https://github.com/dmcable/spacexr) | :heavy_check_mark: | human | [GPL 3.0](https://github.com/dmcable/spacexr/blob/master/LICENSE) | Cable, D.M., Murray, E., Zou, L.S. et al. Robust decomposition of cell type mixtures in spatial transcriptomics. Nat Biotechnol 40, 517–526 (2022). <https://doi.org/10.1038/s41587-021-00830-w>|
| [CARD](https://github.com/YingMa0107/CARD) | :heavy_check_mark: | human | [GPL-3.0](https://github.com/YingMa0107/CARD/blob/master/LICENSE.md) | Ying Ma, and Xiang Zhou (2022). Spatially informed cell type deconvolution for spatial transcriptomics. <https://doi.org/10.1038/s41587-022-01273-7> |
