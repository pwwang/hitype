% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hitype_score.R
\name{hitype_assign_level}
\alias{hitype_assign_level}
\title{Assign cell types based on ScType scores}
\usage{
hitype_assign_level(clusters, scores, threshold, top = 10)
}
\arguments{
\item{clusters}{A named vector of original cluster assignments
(names - cell names, values - cluster assignments)}

\item{scores}{A list of matrices of cell type scores for each level}

\item{threshold}{A threshold for low confidence cell type assignment
The cell types are only assigned for cells with scores higher than the
threshold * ncells.
(0 - 1, default 0.05)}

\item{top}{The number of top cell types to assign for each cluster}
}
\value{
A data from of top cell type assignments
}
\description{
Assign cell types based on ScType scores
}
\keyword{internal}
