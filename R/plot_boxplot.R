#' Compare Mean Results with Boxplots
#'
#' Each colored point is one selected column's mean across spots. Each box
#' summarizes these means within a group, not variation across spots or a
#' confidence interval. Its center line is the median of the column means.
#'
#' @param spe A `SpatialExperiment`.
#' @param column Names of numeric `colData` columns to compare.
#' @param groups Named character vector mapping column names to groups on the
#'   x-axis (e.g. cell types). Order follows first appearance in `column`.
#' @param sources Named character vector mapping column names to sources shown
#'   by point color (e.g. methods). Each group/source pair must be unique.
#' @param title,xlab Optional plot title and x-axis label.
#' @param ylab Y-axis label, default "Mean value".
#' @param colors Optional named character vector mapping sources to colors.
#' @param point_size Size of the points in the plot and legend, in ggplot2 size units.
#' @return A ggplot object. Its `data` contains `column`, `group`, `source`, and
#'   `mean`. Further styling, including percentage axes, can be added with ggplot2.
#' @details Values are used without normalization or aggregation of cell types.
#'   Spots with any nonfinite selected value are omitted from all means, with a
#'   warning. Constant columns are retained. No statistical tests are performed.
#' @export
plot_boxplot <- function(spe, column, groups, sources, title = NULL, xlab = NULL,
                         ylab = "Mean value", colors = NULL, point_size = 3) {
  if (!methods::is(spe, "SpatialExperiment")) stop("spe must be a SpatialExperiment.")
  if (!is.character(column) || !length(column) || anyNA(column) || anyDuplicated(column)) {
    stop("column must contain distinct column names.")
  }
  if (any(!column %in% names(SummarizedExperiment::colData(spe)))) {
    stop("One or more selected columns do not exist in colData(spe).")
  }
  selected <- SummarizedExperiment::colData(spe)[, column, drop = FALSE]
  if (!all(vapply(selected, function(x) is.numeric(x) && is.null(dim(x)), logical(1)))) {
    stop("Selected columns must be numeric vectors.")
  }
  check_mapping <- function(x, keys, argument) {
    if (!is.character(x) || is.null(names(x)) || anyNA(x) ||
        anyNA(names(x)) || anyDuplicated(names(x)) || any(!keys %in% names(x))) {
      stop(argument, " must be a named character vector covering all selected names.")
    }
  }
  check_mapping(groups, column, "groups")
  check_mapping(sources, column, "sources")
  groups <- unname(groups[column])
  sources <- unname(sources[column])
  if (anyDuplicated(data.frame(group = groups, source = sources))) {
    stop("Each group/source pair must correspond to a single column.")
  }
  if (!is.null(colors)) check_mapping(colors, unique(sources), "colors")
  m <- as.matrix(selected)
  valid <- rowSums(!is.finite(m)) == 0L
  if (any(!valid)) warning("Excluded ", sum(!valid), " spots with nonfinite selected values.", call. = FALSE)
  if (!any(valid)) stop("At least one complete spot is required.")
  data <- data.frame(column = column,
    group = factor(groups, levels = unique(groups)),
    source = factor(sources, levels = unique(sources)),
    mean = unname(colMeans(m[valid, , drop = FALSE])))
  group <- source <- mean <- NULL # ggplot data variables
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = group, y = mean)) +
    ggplot2::geom_boxplot(ggplot2::aes(group = group), width = .65,
      outlier.shape = NA, fill = "grey92", color = "grey35", linewidth = .5) +
    ggplot2::geom_point(ggplot2::aes(color = source, group = source),
      position = ggplot2::position_dodge(width = .48), size = point_size) +
    ggplot2::labs(title = title, x = xlab, y = ylab, color = "Source") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())
  if (!is.null(colors)) plot <- plot + ggplot2::scale_color_manual(values = colors)
  plot
}
