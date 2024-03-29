% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train_weights.R
\name{run_weights_on_test_data}
\alias{run_weights_on_test_data}
\title{Run compiled weights on test data}
\usage{
run_weights_on_test_data(db, exprs, clusters, scaled, test_data_idx)
}
\arguments{
\item{db}{The weights data frame}

\item{exprs}{The expression matrix, or a seurat object
(rows: genes, columns: samples/cells)}

\item{clusters}{A named vector of cluster ids
If \code{exprs} is a seurat object, this is ignored. The cluster ids are
taken from the seurat object.}

\item{scaled}{Whether the expression matrix is scaled}

\item{test_data_idx}{The row names of the test data}
}
\description{
Run compiled weights on test data
}
\keyword{internal}
