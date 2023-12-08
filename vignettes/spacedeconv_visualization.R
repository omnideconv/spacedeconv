## ----deconvolute, error=FALSE, message=FALSE, warning=FALSE-------------------
library(spacedeconv)
library(SpatialExperiment)
data("spatial_data_3")

spe <- preprocess(spatial_data_3)
spe <- spacedeconv::normalize(spe, method = "cpm")


deconv <- spacedeconv::deconvolute(
  spe,
  method = "epic",
  assay_sp = "cpm"
)

## ----message=FALSE, warning=FALSE, fig.width=7, fig.height=7------------------
spacedeconv::plot_celltype(deconv, cell_type = "epic_B.cell", density = F, smooth = T)

## ----warning=FALSE, message=FALSE, fig.height=7, fig.width=7------------------
spacedeconv::plot_gene(deconv, gene = "AKT1", density = F, smooth = T)

## ----plotUMIcount, fig.height=7, fig.width=7----------------------------------
spacedeconv::plot_umi_count(deconv, transform_scale = "sqrt", density = F, smooth = T)

## ----plotMostAbundant---------------------------------------------------------
spacedeconv::plot_most_abundant(deconv,
  method = "epic",
  density = F,
  title_size = 25,
  legend_size = 15,
  font_size = 10,
  remove = c("epic_uncharacterized.cell")
)

## ----plotCelltype-------------------------------------------------------------
spacedeconv::plot_celltype_presence(deconv,
  cell_type = "epic_B.cell",
  threshold = 0.05,
  show_image = T
)

## ----plotComparison-----------------------------------------------------------
spacedeconv::plot_comparison(deconv,
  cell_type_1 = "epic_T.cell.CD4.",
  cell_type_2 = "epic_T.cell.CD8.",
  palette = "Vik",
  density = F,
  smooth = T,
)

## ----availablePalettes, fig.width=10------------------------------------------
colorspace::hcl_palettes(plot = TRUE)

## ----offset_rotation, eval = FALSE--------------------------------------------
#  plot_umi_count(deconv, offset_rotation = T) # rotate hexagons

## -----------------------------------------------------------------------------
spacedeconv::plot_umi_count(deconv, transform_scale = "sqrt", density = FALSE)

## -----------------------------------------------------------------------------
spacedeconv::plot_celltype(deconv,
  cell_type = "epic_B.cell",
  smooth = T,
  smoothing_factor = 1.5,
  density = FALSE
)

## -----------------------------------------------------------------------------
spacedeconv::plot_celltype(deconv,
  cell_type = "epic_B.cell",
  smooth = TRUE,
  density = TRUE,
  title_size = 20,
  legend_size = 15,
  font_size = 10
)

## ----eval= FALSE--------------------------------------------------------------
#  spacedeconv::plot_celltype(deconv,
#    cell_type = "epic_B.cell",
#    smooth = TRUE,
#    density = TRUE,
#    save = TRUE,
#    path = "~/project/results",
#    png_width = 1000,
#    png_height = 750
#  )

## ----aggregation, eval=FALSE--------------------------------------------------
#  spe <- aggregate(spe, cell_type_1, cell_type_2, newName)
