# SpaceDeconv

Unified interface for the deconvolution of spatial transcriptomics data

## Installation 

``` r
install.packages("devtools")

# quick install, methods will be installed upon request
devtools::install_github("omnideconv/SpaceDeconv")

# full installation 
devtools::install_github("omnideconv/SpaceDeconv", dependencies = TRUE)
``` 

## Usage 

### Load Spatial Dataset
You can easily load 10xVisium Data by providing the spaceranger output folder. It is further possible to run SpaceDeconv with manually created SpatialExperiments. See the SpatialExperiment [Documentation](https://github.com/drighelli/SpatialExperiment) for further details. 
``` r
spe <- SpatialExperiment::read10xVisium("path_to_spaceranger_output")
``` 

### Normalization 
SpaceDeconv offers an addition function for convenient normalization of SpatialExperiments. The normalization is saved in an additional assay, so make sure the correct data is used during deconvolution by providing the desired assay with the parameters `assay_sc` and `assay_sp`
``` r
spe <- SpaceDeconv::normalize(spe, method="cpm")

```

### Build a Signature Matrix 
``` r
signature <- SpaceDeconv::build_model(single_cell_object, cell_type_col = "cell_ontology_class", method = "spotlight", assay_sc="cpm")
```
### Deconvolution
To perform a deconvolution a SpatialExperiment object is required. Some methods additionally require a cell-type specific reference signature which can be calculated by `SpaceDeconv::build_model()`. By default the deconvolution results are added to the SpatialExperiment object for convenient visualization. You can obtain the results in table form by setting `return_object=FALSE`.
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
