

print_info <- function(sce = NULL, spe = NULL) {
  # check for correct class

  # general info
  width <- unname(unlist(options("width")))

  message(rep("=", floor((width - 10) / 2)), " spacedeconv ", rep("=", floor((width - 10) / 2)))


  # exclamation, rocket, check, cross mark, ok, hourglass, entry
  # sce
  if (!is.null(sce)) {
    message()
    message("=== Single Cell ===")
    message("Assays: ", paste0(assayNames(sce), collapse = ", "))
    message("Genes: ", nrow(sce))

    genes0 <- sum(rowSums(counts(sce)) == 0)
    percentGenes0 <- round(genes0 / nrow(sce) * 100, 2)
    message("> without expression: ", genes0, " (", percentGenes0, "%)")

    message("Cells: ", ncol(sce))
    umi <- colSums(counts(sce))
    cells0 <- sum(umi == 0)
    percentCells0 <- round(cells0 / ncol(sce) * 100, 2)
    message("> without expression: ", cells0, " (", percentCells0, "%)")
    message("Umi count range: ", min(umi), "-", max(umi))


    if (is.null(rownames(sce))) {
      message(emo::ji("cross mark"), " no rownames set")
    } else {
      message(emo::ji("check"), " rownames set")
    }

    if (is.null(colnames(sce))) {
      message(emo::ji("cross mark"), " no colnames set")
    } else {
      message(emo::ji("check"), " colnames set")
    }
  }


  # spe

  if (!is.null(spe)) {
    message()
    message("=== Spatial ===")
    message("Assays: ", paste0(assayNames(spe), collapse = ", "))
    message("Genes: ", nrow(spe))

    genes0 <- sum(rowSums(counts(spe)) == 0)
    percentGenes0 <- round(genes0 / nrow(spe) * 100, 2)
    message("> without expression: ", genes0, " (", percentGenes0, "%)")

    message("Spots: ", ncol(spe))
    umi <- colSums(counts(spe))
    cells0 <- sum(umi == 0)
    percentCells0 <- round(cells0 / ncol(spe) * 100, 2)
    message("> without expression: ", cells0, " (", percentCells0, "%)")
    message("Umi count range: ", min(umi), "-", max(umi))


    if (is.null(rownames(spe))) {
      message(emo::ji("cross mark"), " no rownames set")
    } else {
      message(emo::ji("check"), " rownames set")
    }

    if (is.null(colnames(spe))) {
      message(emo::ji("cross mark"), " no colnames set")
    } else {
      message(emo::ji("check"), " colnames set")
    }
  }
}

print_info(sce = sce, spe = spe)
