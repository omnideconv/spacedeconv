# SpaceDeconv

Unified interface for the deconvolution of spatial transcriptomics data

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
signature <- SpaceDeconv::build_model(single_cell_object, cell_type_col = "cell_ontology_class", method = "spotlight", assay_sc="cpm")
```
### Deconvolution
To perform a deconvolution a SpatialExperiment object is required. Some methods additionally require a cell-type specific reference signature which can be calculated by `SpaceDeconv::build_model()`. By default the deconvolution results are added to the SpatialExperiment object to simplify the visualization. You can obtain the results in table form by setting `return_object=FALSE`.
```r
# save the results to an annotated SpatialExperiment
result <- SpaceDeconv::deconvolute(spatial_object, signature, method = "spotlight")

# return deconvolution results as a table
result <- SpaceDeconv::deconvolute(spatial_object, signature, method = "spotlight", return_object = FALSE)
```


### Visualization 
SpaceDeconv includes multiple visualization functions. 
```r
# sample does refer to the first column of ColData(spe)
# for cell_type input a celltype present in the deconvolution result
plot_celltype(spe, sample="sammple01", cell_type"B.cells")

# threshold changes the minimum cell type fraction for a cell to be considered present in a specific spot
plot_cells_per_spot(spe, threshold=0.01)
```
