library("Orthology.eg.db")
library("org.Mm.eg.db")
library("org.Hs.eg.db")

#' map human to mouse genes
#'
#' @param humangenes a vector containing human genes to be converted to mouse
#'
#' @export
#'

convert_human_to_mouse <- function(humangenes) {
  # Function to process individual or complex gene names
  processGene <- function(gene) {
    # Split gene names if they are in complex form (e.g., "GENE1_GENE2")
    splitGenes <- unlist(strsplit(gene, "_"))

    # Process each gene symbol individually
    mappedSymbols <- sapply(splitGenes, function(x) {
      gns <- tryCatch(
        {
          mapIds(org.Hs.eg.db, x, "ENTREZID", "SYMBOL")
        },
        error = function(e) NA
      )
      if (!is.na(gns)) {
        mapped <- tryCatch(
          {
            select(Orthology.eg.db, gns, "Mus_musculus", "Homo_sapiens")
          },
          error = function(e) {
            return(NA)
          }
        )
        if (!is.na(mapped$Mus_musculus)) {
          msymb <- tryCatch(
            {
              mapIds(org.Mm.eg.db, as.character(mapped$Mus_musculus), "SYMBOL", "ENTREZID")
            },
            error = function(e) NA
          )
          return(msymb)
        }
      }
      return(NA)
    })
    # Combine the mapped symbols with underscores if they were a complex
    paste(mappedSymbols, collapse = "_")
  }

  suppressWarnings({
    mappedGenes <- sapply(humangenes, processGene)
  })

  out <- data.frame(Human_symbol = humangenes, Mouse_symbol = mappedGenes, stringsAsFactors = FALSE)
  return(out)
}


#' map mouse to human genes
#'
#' @param mousegenes a vector containing mouse genes to be converted to human
#'
#' @export
#'

convert_mouse_to_human <- function(mousegenes) {
  # Function to process individual or complex gene names
  processGene <- function(gene) {
    # Split gene names if they are in complex form (e.g., "GENE1_GENE2")
    splitGenes <- unlist(strsplit(gene, "_"))

    # Process each gene symbol individually
    mappedSymbols <- sapply(splitGenes, function(x) {
      gns <- tryCatch(
        {
          mapIds(org.Mm.eg.db, x, "ENTREZID", "SYMBOL")
        },
        error = function(e) NA
      )
      if (!is.na(gns)) {
        mapped <- tryCatch(
          {
            select(Orthology.eg.db, gns, "Homo_sapiens", "Mus_musculus")
          },
          error = function(e) {
            return(NA)
          }
        )
        if (!is.na(mapped$Homo_sapiens)) {
          hsymb <- tryCatch(
            {
              mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens), "SYMBOL", "ENTREZID")
            },
            error = function(e) NA
          )
          return(hsymb)
        }
      }
      return(NA)
    })
    # Combine the mapped symbols with underscores if they were a complex
    paste(mappedSymbols, collapse = "_")
  }

  suppressWarnings({
    mappedGenes <- sapply(mousegenes, processGene)
  })
  out <- data.frame(Mouse_symbol = mousegenes, Human_symbol = mappedGenes, stringsAsFactors = FALSE)
  return(out)
}
