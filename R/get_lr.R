#' Get Ligand Expression
#'
#' @description This function computes the expression values of a specified ligand (gene or gene complex) from a given expression dataset.
#' @param gene_pair A string representing a gene or a gene complex. Gene complexes are indicated by genes separated by underscores.
#' @param cpm_df A dataframe containing gene expression data, where row names are gene names and columns represent different samples or spots.
#' @return A numeric vector containing the expression values of the specified ligand across the samples or spots in `cpm_df`. If the ligand is a gene complex, it returns the mean expression value of the component genes. If any component gene is missing or has zero expression, it returns a vector of zeros.
get_ligand_expression <- function(gene_pair, cpm_df) {
  first_gene <- sub("~.*", "", gene_pair) # Extract the first gene/gene complex before "~"

  # check if its a complex or not
  if (grepl("_", first_gene)) {
    sub_genes <- unlist(strsplit(first_gene, "_"))
    # Check if all sub_genes exist in rownames(cpm_df) and have non-zero expression values
    if (all(sub_genes %in% rownames(cpm_df)) &&
      all(cpm_df[sub_genes, ] != 0)) {
      # Calculate the mean expression values for each sub-gene across spots
      expression_values <- apply(cpm_df[sub_genes, ], 2, mean)
    } else {
      # If any sub_gene does not exist or has zero expression, return a vector of zeros
      expression_values <- rep(0, ncol(cpm_df))
    }
  } else {
    expression_values <- cpm_df[first_gene, ]
    # return(expression_values)
  }

  expression_values
}


#' Get Receptor Expression
#'
#' @description This function calculates the expression values of a specified receptor (gene or gene complex) from a given expression dataset.
#' @param gene_pair A string representing a gene or a gene complex, where gene complexes are indicated by genes separated by underscores, and the receptor gene is specified after a "~".
#' @param cpm_df A dataframe containing gene expression data, with row names as gene names and columns as different samples or spots.
#' @return A numeric vector of the expression values of the specified receptor across the samples or spots in `cpm_df`. If the receptor is a gene complex, the function returns the mean expression value of the component genes. If any component gene is missing or has zero expression, a vector of zeros is returned.
get_receptor_expression <- function(gene_pair, cpm_df) {
  second_gene <- sub(".*~", "", gene_pair) # Extract the second gene/gene complex after "~"

  # check if its a complex or not
  if (grepl("_", second_gene)) {
    sub_genes <- unlist(strsplit(second_gene, "_"))
    # Check if all sub_genes exist in rownames(cpm_df) and have non-zero expression values
    if (all(sub_genes %in% rownames(cpm_df)) &&
      all(cpm_df[sub_genes, ] != 0)) {
      # Calculate the mean expression values for each sub-gene across spots
      expression_values <- apply(cpm_df[sub_genes, ], 2, mean)
    } else {
      # If any sub_gene does not exist or has zero expression, return a vector of zeros
      expression_values <- rep(0, ncol(cpm_df))
    }
  } else {
    expression_values <- cpm_df[second_gene, ]
    # return(expression_values)
  }

  expression_values
}

#' Compute L-R score for each spot
#'
#' @param spe SpatialExperiment object
#' @param resource "consensus" table from Omnipath as default, the user can provide a data frame containing L-R pairs. The data frame should contain at least the following two columns: source_genesymbol = ligands, target_genesymbol = receptors
#' @param method mathematical approach to compute L-R scores. Options: min, and product L-R, default is product
#' @param organism choose the organism to be considered, default human, options: human or mouse
#'
#' @export
get_lr <- function(spe,
                   resource = "Consensus",
                   method = "min",
                   organism = "human") {
  startTime <- Sys.time()


  # Checking Parameter
  if (is.null(spe)) {
    stop("Please provide a SpatialExperiment object")
  }

  if (!is.null(method) && !(method %in% c("min", "product"))) {
    stop("The method provided is not supported. Please provide one of the following methods to calculate the L-R scores: min, product")
  }

  if (!is.null(organism) && !(organism %in% c("human", "mouse"))) {
    stop("The organism provided is not supported. Please provide one of the following organisms: human, mouse")
  }

  # Check ENSEMBL identifiers
  if (any(grepl("^ENS[A-Z]+[0-9]+", rownames(SummarizedExperiment::assay(object, assay))))) {
    cli::cli_alert_warning("Warning: ENSEMBL identifiers detected in gene names")
    stop("Please convert ENSEMBL to Gene ID")
  }

  # normalize the data and extraxt the cpm values
  # spe <- spacedeconv::preprocess(spe)
  # spe <- spacedeconv::normalize(spe, method = "cpm")
  if (!"cpm" %in% assayNames(spe)) {
    stop("Please provide an object with cpm normalized values")
  }

  # extract the cpm matrix
  # first as matrix and then as df
  cpm_matrix <- as.matrix(spe@assays@data$cpm)
  cpm_df <- as.data.frame(cpm_matrix)
  rm(cpm_matrix)

  # check if mouse organism is provided
  if (organism == "mouse") {
    # if resource not defined or resource defined as "Consensus" then take the LIANA Consensus
    if (is.null(resource) || resource == "Consensus") {
      # import the omnipath intercell network
      resource <- import_intercell_network()
      cis <- c(
        "cell_surface_ligand",
        "ecm",
        "ligand",
        ""
      )
      cit <- c(
        "receptor",
        "adhesion",
        ""
      )
      sc <- c(
        "CellPhoneDB",
        "CellChatDB",
        "ICELLNET",
        "connectomeDB2020",
        "CellTalkDB"
      )
      # perform filtering based on the arguments above
      resource <- resource[resource$category_intercell_source %in% cis, ]
      resource <- resource[resource$category_intercell_target %in% cit, ]
      resource <- resource[grep(paste(sc, collapse = "|"), resource$sources), ]
      resource <- resource[resource$omnipath, ]
      resource <- resource[, c("source_genesymbol", "target_genesymbol")]
      # create pairs
      resource$pairs <- paste(resource$source_genesymbol, resource$target_genesymbol, sep = "~")
      # remove duplicates
      resource <- resource[!duplicated(resource$pairs), ]
      # conversion from human to mouse
      # ligands
      cli::cli_progress_step(msg = "Converting the human ligands to mouse", msg_done = "Conversion completed!")
      suppressWarnings({
        mligands <- convert_human_to_mouse(humangenes = resource$source_genesymbol)
      })
      cli::cli_progress_done()

      # receptors
      cli::cli_progress_step(msg = "Converting the human receptors to mouse", msg_done = "Conversion completed!")
      suppressWarnings({
        mreceptors <- convert_human_to_mouse(humangenes = resource$target_genesymbol)
      })
      cli::cli_progress_done()

      # Create a temporary data frame with both Mouse_symbol columns
      temp_df <- data.frame(mligands_M = mligands$Mouse_symbol, mreceptors_M = mreceptors$Mouse_symbol)

      # Create a logical vector indicating rows without "NA" in either column
      valid_rows <- !(temp_df$mligands_M == "NA" | temp_df$mreceptors_M == "NA")

      # Use the logical vector to subset the data
      temp_df <- temp_df[valid_rows, ]

      # Create the pairs by pasting the columns together
      temp_df$pairs <- paste(temp_df$mligands_M, temp_df$mreceptors_M, sep = "~")

      resource <- temp_df
    } else {
      if (!is.data.frame(resource)) {
        stop("The 'resource' parameter must be a data frame.")
      } else {
        # Check if both columns exist in the resource dataframe
        if ("source_genesymbol" %in% colnames(resource) & "target_genesymbol" %in% colnames(resource)) {
          # Both columns exist
          resource$pairs <- paste(resource$source_genesymbol, resource$target_genesymbol, sep = "~")
        } else {
          # At least one of the columns doesn't exist
          stop("At least one of the columns does not exist in the data frame provided or has a different name. The column names should be source_symbol for the ligands and target_symbol for the receptors")
        }
      }
    }
  } else {
    ######### HUMAN #############

    # check if the organism is human
    if (is.null(organism) || organism == "human") {
      # if nothing is provided then take the LIANA consensus
      if (is.null(resource) || resource == "Consensus") {
        cli::cli_alert_info("Using the Consensus reference")
        resource <- import_intercell_network()
        cis <- c(
          "cell_surface_ligand",
          "ecm",
          "ligand",
          ""
        )
        cit <- c(
          "receptor",
          "adhesion",
          ""
        )
        sc <- c(
          "CellPhoneDB",
          "CellChatDB",
          "ICELLNET",
          "connectomeDB2020",
          "CellTalkDB"
        )
        # perform filtering based on the arguments above
        resource <- resource[resource$category_intercell_source %in% cis, ]
        resource <- resource[resource$category_intercell_target %in% cit, ]
        resource <- resource[grep(paste(sc, collapse = "|"), resource$sources), ]
        resource <- resource[resource$omnipath, ]
        # Select only the columns 'source_genesymbol' and 'target_genesymbol'
        resource <- resource[, c("source_genesymbol", "target_genesymbol")]
        # create pairs
        resource$pairs <- paste(resource$source_genesymbol, resource$target_genesymbol, sep = "~")
        # remove duplicates
        resource <- resource[!duplicated(resource$pairs), ]
      } else {
        if (!is.data.frame(resource)) {
          stop("The 'resource' parameter must be a data frame.")
        } else {
          # Check if both columns exist in the resource dataframe
          if ("source_genesymbol" %in% colnames(resource) & "target_genesymbol" %in% colnames(resource)) {
            # Both columns exist
            resource$pairs <- paste(resource$source_genesymbol, resource$target_genesymbol, sep = "~")
          } else {
            # At least one of the columns doesn't exist
            stop("At least one of the columns does not exist in the data frame provided or has a different name. The column names should be source_symbol for the ligands and target_symbol for the receptors")
          }
        }
      }
    } else {
      stop("The provided organism is not supported. Choose between human or mouse")
    }
  }

  cli::cli_progress_step("Extracting Expression Data")

  # prepare the ligand and receptor data frames which will contain the cpm values for each
  # ligands
  cpm_ligands <- data.frame(matrix(NA, nrow = nrow(resource), ncol = ncol(cpm_df)))

  # Set row names based on the "pairs" column
  row.names(cpm_ligands) <- resource$pairs
  colnames(cpm_ligands) <- colnames(cpm_df)
  # the same for the receptors
  cpm_receptors <- cpm_ligands

  # fill the ligand table with the cpm values
  # Function to get expression values for the the ligands (first elementa of the L-R pairing)
  cli::cli_progress_step("Extracting Ligand Expression")

  # Fill in cpm_ligands directly
  for (row_name in rownames(cpm_ligands)) {
    gene_pair <- row_name
    expression_values <- get_ligand_expression(gene_pair, cpm_df)
    cpm_ligands[row_name, ] <- as.numeric(expression_values)
  }

  # receptors
  cli::cli_progress_step("Extracting Receptor Expression")

  # Fill in cpm_receptors directly
  for (row_name in rownames(cpm_receptors)) {
    gene_pair <- row_name
    expression_values <- get_receptor_expression(gene_pair, cpm_df)
    cpm_receptors[row_name, ] <- as.numeric(expression_values)
  }

  cli::cli_alert_info("Removing NA values from both tables")
  # filtering
  # clean the ligands table
  cpm_ligands_clean <- na.omit(cpm_ligands)

  cli::cli_progress_step("Preparing Expression Data for computing L-R score")

  cpm_receptors_clean <- na.omit(cpm_receptors)

  # find the common rows
  common_rows <- intersect(rownames(cpm_ligands_clean), rownames(cpm_receptors_clean))

  # Subset the data frames to only include common rows
  cpm_ligands_subset <- cpm_ligands_clean[common_rows, ]
  cpm_receptors_subset <- cpm_receptors_clean[common_rows, ]

  cli::cli_progress_step(msg = "Calculating L-R score")
  # choose the mathematical operation
  # if no method provide then the default (min) should be used
  if (is.null(method) || method == "min") {
    # calculate the minimum expression value per spot as default
    result <- pmin(cpm_ligands_subset, cpm_receptors_subset)
  } else if (method == "product") {
    # Perform element-wise multiplication
    result <- cpm_ligands_subset * cpm_receptors_subset
  } else {
    stop("The method provided is not supported. Please provide one of the following: product, min")
  }

  result <- t(result)

  # replace the colnames and add a lr prefix
  colnames(result) <- make.names(colnames(result))
  colnames(result) <- paste("lr_", colnames(result), sep = "")

  # add the result as additional columns to the colData of the spe object
  colData(spe) <- cbind(colData(spe), result)

  cli::cli_progress_done()

  endTime <- Sys.time()

  # Print the execution time more clearly
  cli::cli_alert_success(paste("Finished computation in", format(difftime(endTime, startTime, units = "mins"), digits = 4), " minutes"))

  return(spe)
}
