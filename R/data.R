# Built-in gene sets and marker databases shipped with hitype
#
# These objects live in data/ and are documented here so that R CMD check
# does not flag them as undocumented datasets.

#' Gene sets from ScType (short version)
#'
#' The short marker database from
#' \href{https://github.com/IanevskiAleksandr/sc-type}{ScType}, bundled for
#' convenience. It can be passed directly to \code{\link{gs_prepare}}:
#'
#' \code{gs <- gs_prepare(hitypedb_short, tissue_type = "Immune system")}
#'
#' @format A data frame with columns commonly used by \code{\link{gs_prepare}}:
#'   \code{cellName}, \code{tissueType}, \code{geneSymbolmore1},
#'   \code{geneSymbolmore2} and \code{level}.
#' @source \url{https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx}
"hitypedb_short"

#' Gene sets from ScType (full version)
#'
#' The full marker database from
#' \href{https://github.com/IanevskiAleksandr/sc-type}{ScType}, bundled for
#' convenience. It can be passed directly to \code{\link{gs_prepare}}, e.g.
#' \code{gs_prepare(hitypedb_full, tissue_type = "Immune system")}.
#'
#' @format A data frame with columns commonly used by \code{\link{gs_prepare}}:
#'   \code{cellName}, \code{tissueType}, \code{geneSymbolmore1},
#'   \code{geneSymbolmore2} and \code{level}.
#' @source \url{https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx}
"hitypedb_full"

#' Gene sets with weights trained from the PBMC 3k dataset
#'
#' Marker gene sets for the 9 cell types of the Seurat PBMC 3k tutorial
#' (\code{Naive CD4+ T}, \code{CD14+ Mono}, \code{Memory CD4+}, \code{B},
#' \code{CD8+ T}, \code{FCGR3A+ Mono}, \code{NK}, \code{DC} and
#' \code{Platelet}), with weights trained on the dataset itself. It is
#' intended for demonstration: pass it to \code{\link{gs_prepare}} like any
#' other marker data frame (e.g. \code{gs_prepare(hitypedb_pbmc3k)}).
#'
#' @format A data frame in the marker database format used by
#'   \code{\link{gs_prepare}}: \code{cellName}, \code{geneSymbolmore1},
#'   \code{geneSymbolmore2} and \code{level}, where \code{geneSymbolmore1}
#'   carries \code{+}/\code{-} suffixes encoding the trained weights.
#' @source Trained with \code{\link{train_weights}} on the PBMC 3k dataset
#'   (Seurat guided-clustering tutorial, resolution 0.5).
"hitypedb_pbmc3k"