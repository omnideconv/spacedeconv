% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{get_palette_type}
\alias{get_palette_type}
\title{Determine the Type of a Given Color Palette}
\usage{
get_palette_type(palette)
}
\arguments{
\item{palette}{A character string specifying the name of the color palette.
The name should be a valid palette name from either the
\code{colorspace} or \code{RColorBrewer} package.}
}
\value{
A character string representing the type of the palette.
Possible return values are "sequential", "diverging", or "qualitative".
The function will stop and display an error message if an unknown palette
name is provided.
}
\description{
This function identifies the type of a color palette, supporting palettes
from both the \code{colorspace} and \code{RColorBrewer} packages. It is used to
streamline the process of palette type identification for use in plotting functions.
}
