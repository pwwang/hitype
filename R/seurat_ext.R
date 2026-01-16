#' Run hitype_assign for a Seurat object
#'
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat Idents
#'
#' @rdname RunHitype
#'
#' @param object Seurat object
#' @param gs The gene list prepared by \code{\link{gs_prepare}}
#' @param fallback A fallback cell type if no cell type is assigned
#' @param threshold A threshold for low confidence cell type assignment
#'  The cell types are only assigned for cells with scores higher than the
#'  threshold.
#'  (0 - 1, default 0.05)
#' @param level_weights The weights for each level of the hierarchy to calculate
#'  the final cell type score
#'  It should be either a numeric vector of length equal to the number of levels
#'  or a single numeric value to be used for all levels
#'  It can also be a function that takes the levels as input and returns a
#'  numeric vectors as the weights.
#' @param make_unique Whether to make the cell type names unique
#' @param layer The layer to use for `GetAssayData`
#' @param assay The assay to use for `GetAssayData`
#' @param scaled Whether the data from `GetAssayData` is scaled
#' @return The Seurat object with the cell types (named `hitype`) added to the
#'  metadata
#' @export
RunHitype <- function(object, ...) {
    UseMethod(generic = "RunHitype", object = object)
}

#' @rdname RunHitype
#' @export
RunHitype.default <- function(object, ...) {
    stop("RunHitype is not implemented for this object type")
}

#' @rdname RunHitype
#' @export
RunHitype.Seurat <- function(
    object,
    gs = NULL,
    fallback = "Unknown",
    threshold = 0.05,
    level_weights = function(l) 1 / (10 ^ (l - 1)),
    make_unique = FALSE,
    layer = "data",
    assay = NULL,
    scaled = FALSE,
    ...
) {
    scores <- hitype_score(
        Seurat::GetAssayData(object, layer = layer, assay = assay),
        gs = gs,
        scaled = scaled
    )
    clusters <- Seurat::Idents(object)
    cell_types <- hitype_assign(
        clusters,
        scores = scores,
        gs = gs,
        fallback = fallback,
        threshold = threshold
    )
    # Level, Cluster, CellType, Score
    cell_types <- summary(
        cell_types,
        level_weights = level_weights,
        make_unique = make_unique
    )
    # Add to metadata
    object@meta.data$hitype <- cell_types[
        match(clusters, cell_types$Cluster),
        "CellType",
        drop = TRUE
    ]
    object
}
