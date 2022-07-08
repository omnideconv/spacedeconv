# SpaceDeconv

Unified interface for the deconvolution of spatial transcriptomics data

## Installation 

``` r
install.packages("devtools")

# quick install, methods will be installed upon request
devtools::install_github("omnideconv/SpaceDeconv")
``` 

## Usage 

### 1. Build a Signature Matrix 

``` r
signature <- SpaceDeconv::build_model(single_cell_object, cell_type_col = "cell_ontology_class", method = "spotlight")
```

### 2. Deconvolute 

```r
result <- SpaceDeconv::deconvolute(spatial_object, signature, method = "spotlight")
```
