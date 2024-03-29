% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hitype_score.R
\name{hitype_score_level}
\alias{hitype_score_level}
\title{Calculate ScType scores and assign cell types for one level}
\usage{
hitype_score_level(z, gs_level)
}
\arguments{
\item{z}{Z-scaled expression matrix
(rownames - gene names, colnames - cell names)}

\item{gs_level}{One level of gene sets prepared by \code{\link{gs_prepare}}}
}
\value{
A matrix of cell type scores for each cell
(rownames - cell types, colnames - cell names)
}
\description{
Calculate ScType scores and assign cell types for one level
}
\keyword{internal}
