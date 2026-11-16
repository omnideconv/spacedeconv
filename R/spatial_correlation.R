#' Create a Correlation Heatmap
#'
#' Correlates numeric columns of `colData(spe)` across spots and draws a full
#' heatmap in the supplied column order, without clustering.
#'
#' @param spe A `SpatialExperiment`.
#' @param column Character vector of at least two numeric `colData` column names,
#'   in the desired display order.
#' @param cor_method Correlation method: "pearson", "spearman", or "kendall".
#' @param categories Optional named character vector mapping column names to
#'   categories, shown as color strips above and beside the heatmap.
#' @param category_colors Optional named vector mapping categories to colors.
#'   When NULL, colors are assigned automatically in category appearance order.
#' @param labels Optional named character vector mapping column names to row
#'   labels. Defaults to the column names.
#' @param title Optional heatmap title; NULL omits it.
#' @param label_title Heading above the row labels, defaulting to
#'   "Source | Feature | Mean"; NULL omits it.
#'   For example, "Method | Cell type | Mean relative abundance". This is free
#'   text and does not affect the values or labels.
#' @param annotation_legend_param Named list of category legend settings passed
#'   to [ComplexHeatmap::HeatmapAnnotation()], such as `nrow` or `labels_gp`.
#' @param label_title_gp Graphical parameters for the heading above row labels.
#' @param ... Additional arguments to [ComplexHeatmap::Heatmap()], such as
#'   `row_names_gp`, `column_title_gp`, or `heatmap_legend_param`. These override
#'   styling defaults, but cannot override the matrix, ordering, labels, title,
#'   category annotations, or enable clustering or splitting.
#'
#' @details Values are used as supplied: no normalization, cell-type aggregation,
#'   or percentage conversion is performed. Spots with any nonfinite selected
#'   value are excluded jointly, then constant columns are excluded, with warnings.
#'   At least two spots and two variable columns must remain. No p-values are
#'   calculated. This interface replaces the former corrplot-based interface;
#'   select explicit columns with `column` rather than a method prefix.
#'
#' @return Invisibly, a list containing `corr` (correlation matrix), `heatmap`
#'   (the drawn ComplexHeatmap object), `used_spots`, `excluded_spots` (spot names,
#'   or positional identifiers if unnamed), and `excluded_columns` (constant
#'   column names). The heatmap can be redrawn with `ComplexHeatmap::draw()`.
#' @export
spatialcorr <- function(spe, column, cor_method = "pearson", categories = NULL,
                        category_colors = NULL, labels = NULL, title = NULL,
                        label_title = "Source | Feature | Mean",
                        annotation_legend_param = list(),
                        label_title_gp = grid::gpar(fontsize = 12, fontface = "bold"), ...) {
  if (!methods::is(spe, "SpatialExperiment")) stop("spe must be a SpatialExperiment.")
  if (!is.character(column) || length(column) < 2L || anyNA(column) || anyDuplicated(column)) {
    stop("column must contain at least two distinct column names.")
  }
  if (any(!column %in% names(SummarizedExperiment::colData(spe)))) {
    stop("One or more selected columns do not exist in colData(spe).")
  }
  selected <- SummarizedExperiment::colData(spe)[, column, drop = FALSE]
  if (!all(vapply(selected, function(x) is.numeric(x) && is.null(dim(x)), logical(1)))) {
    stop("Selected columns must be numeric vectors.")
  }
  cor_method <- match.arg(cor_method, c("pearson", "spearman", "kendall"))
  # Match metadata by names, never by the order of the supplied vectors.
  check_mapping <- function(x, keys, argument) {
    if (!is.character(x) || is.null(names(x)) || anyNA(x) ||
        anyNA(names(x)) || anyDuplicated(names(x)) || any(!keys %in% names(x))) {
      stop(argument, " must be a named character vector covering all selected names.")
    }
  }
  if (!is.null(labels)) check_mapping(labels, column, "labels")
  if (!is.null(categories)) check_mapping(categories, column, "categories")
  if (!is.null(category_colors)) {
    if (is.null(categories)) stop("category_colors requires categories.")
    check_mapping(category_colors, unique(categories[column]), "category_colors")
  }
  if (!is.null(label_title) && (!is.character(label_title) || length(label_title) != 1L || is.na(label_title))) {
    stop("label_title must be a single string or NULL.")
  }
  if (!is.list(annotation_legend_param) ||
      (length(annotation_legend_param) && (is.null(names(annotation_legend_param)) ||
       anyNA(names(annotation_legend_param)) || any(names(annotation_legend_param) == "") ||
       anyDuplicated(names(annotation_legend_param))))) {
    stop("annotation_legend_param must be a uniquely named list.")
  }
  if (!inherits(label_title_gp, "gpar")) stop("label_title_gp must be a grid::gpar object.")
  extra <- list(...)
  reserved <- c("matrix", "name", "cluster_rows", "cluster_columns", "row_order",
                "column_order", "row_split", "column_split", "row_km", "column_km",
                "row_labels", "column_labels", "column_title", "top_annotation",
                "left_annotation", "right_annotation")
  if (length(extra) && (is.null(names(extra)) || any(names(extra) == "") || anyDuplicated(names(extra)))) {
    stop("Additional heatmap arguments must have unique names.")
  }
  if (any(names(extra) %in% reserved)) {
    stop("Cannot override fixed heatmap arguments: ", paste(intersect(names(extra), reserved), collapse = ", "))
  }
  m <- as.matrix(selected)
  spots <- colnames(spe)
  if (is.null(spots)) spots <- as.character(seq_len(ncol(spe)))
  valid <- rowSums(!is.finite(m)) == 0L
  if (any(!valid)) warning("Excluded ", sum(!valid), " spots with nonfinite selected values.", call. = FALSE)
  m <- m[valid, , drop = FALSE]
  if (nrow(m) < 2L) stop("At least two complete spots are required.")
  variable <- vapply(seq_len(ncol(m)), function(i) length(unique(m[, i])) > 1L, logical(1))
  excluded_columns <- column[!variable]
  if (length(excluded_columns)) {
    warning("Excluded constant columns: ", paste(excluded_columns, collapse = ", "), call. = FALSE)
  }
  m <- m[, variable, drop = FALSE]
  if (ncol(m) < 2L) stop("At least two nonconstant columns are required.")
  column <- colnames(m)
  r <- stats::cor(m, method = cor_method)
  if (any(!is.finite(r))) stop("Correlation calculation produced nonfinite values.")
  top <- left <- NULL
  if (!is.null(categories)) {
    groups <- unique(unname(categories[column]))
    if (is.null(category_colors)) {
      category_colors <- setNames(grDevices::hcl.colors(length(groups), "Dark 3"), groups)
    }
    ann <- data.frame(Category = unname(categories[column]))
    colors <- list(Category = category_colors[groups])
    legend_params <- list(title_gp = grid::gpar(fontsize = 18, fontface = "bold"),
                          labels_gp = grid::gpar(fontsize = 12))
    legend_params[names(annotation_legend_param)] <- annotation_legend_param
    top <- ComplexHeatmap::HeatmapAnnotation(df = ann, col = colors,
      show_annotation_name = FALSE,
      annotation_legend_param = list(Category = legend_params))
    left <- ComplexHeatmap::rowAnnotation(df = ann, col = colors, show_legend = FALSE,
      show_annotation_name = FALSE)
  }
  legend_name <- switch(cor_method, pearson = "Pearson r", spearman = "Spearman rho", kendall = "Kendall tau")
  args <- list(matrix = r, name = legend_name,
    col = circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "#FFFFFF", "#B2182B")),
    cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = top, left_annotation = left,
    row_labels = if (is.null(labels)) column else unname(labels[column]),
    row_names_max_width = grid::unit(13, "cm"), row_names_gp = grid::gpar(fontsize = 10.5),
    show_row_names = TRUE, show_column_names = FALSE, column_title = if (is.null(title)) character(0) else title,
    column_title_gp = grid::gpar(fontsize = 22), border = TRUE,
    heatmap_legend_param = list(at = c(-1, -.5, 0, .5, 1),
      title_gp = grid::gpar(fontsize = 18, fontface = "bold"), labels_gp = grid::gpar(fontsize = 12)))
  args[names(extra)] <- extra
  if (!is.null(label_title) && isTRUE(args$show_row_names)) {
    # A text annotation retains its heading when the returned heatmap is redrawn.
    heading_gp <- label_title_gp
    label_width <- grid::unit.pmax(
      ComplexHeatmap::max_text_width(args$row_labels, gp = args$row_names_gp),
      ComplexHeatmap::max_text_width(label_title, gp = heading_gp))
    args$show_row_names <- FALSE
    args$right_annotation <- ComplexHeatmap::rowAnnotation(
      Labels = ComplexHeatmap::anno_text(args$row_labels, gp = args$row_names_gp,
        just = "left", location = grid::unit(0, "npc"), width = label_width, show_name = TRUE),
      annotation_label = label_title, annotation_name_side = "top",
      annotation_name_rot = 0, annotation_name_gp = heading_gp)
    args$right_annotation@anno_list$Labels@name_param$x <- grid::unit(0, "npc")
    args$right_annotation@anno_list$Labels@name_param$just <- c("left", "bottom")
  }
  heatmap <- ComplexHeatmap::draw(do.call(ComplexHeatmap::Heatmap, args),
    merge_legends = TRUE, padding = grid::unit(c(2, 2, 2, 10), "mm"))
  invisible(list(corr = r, heatmap = heatmap, used_spots = spots[valid],
                 excluded_spots = spots[!valid], excluded_columns = excluded_columns))
}
