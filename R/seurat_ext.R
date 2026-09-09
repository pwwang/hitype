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
#' @param threshold Confidence threshold as top1/top2 score ratio,
#'  passed to [hitype_assign()]. `NULL` (default) means no filtering.
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
#' @param ident The identity column to use majority voting to assign cell types to clusters
#'  If NULL, return cell-level assignments for each cell.
#'  If "ident", return cluster(identity)-level assignments for each cluster.
#'  if a character, return cluster-level assignments for each cluster based on the specified column in the metadata.
#' @param ... Additional arguments passed to the specific method.
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
    threshold = NULL,
    level_weights = function(l) 1 / (10 ^ (l - 1)),
    make_unique = FALSE,
    layer = "data",
    assay = NULL,
    scaled = FALSE,
    ident = NULL,
    ...
) {
    scores <- hitype_score(
        Seurat::GetAssayData(object, layer = layer, assay = assay),
        gs = gs,
        scaled = scaled
    )
    if (is.null(ident)) {
        # cell-level assignments: give each cell its own "cluster" so that
        # every cell gets its own top cell type
        clusters <- Seurat::Cells(object)
        names(clusters) <- clusters
    } else if (identical(ident, "ident")) {
        clusters <- Seurat::Idents(object)
    } else if (is.character(ident) && length(ident) == 1 && ident %in% colnames(object@meta.data)) {
        # hitype_assign() maps score columns to cells by the cluster names,
        # so the vector must be named (as the Idents branch above is)
        clusters <- object@meta.data[[ident]]
        names(clusters) <- colnames(object)
    } else {
        stop("Invalid ident argument. It should be NULL, 'ident', or a character string corresponding to a column in the metadata.")
    }

    cell_types <- hitype_assign(
        clusters,
        scores = scores,
        gs = gs,
        fallback = fallback,
        threshold = threshold
    )
    # Level, Cluster, CellType, Score
    # make_unique dedupes repeated type names across clusters; it is not
    # applicable to cell-level assignments where each row is one cell
    cell_types <- summary(
        cell_types,
        level_weights = level_weights,
        make_unique = make_unique && !is.null(ident)
    )
    # Add to metadata
    object@meta.data$hitype <- cell_types[
        match(clusters, cell_types$Cluster),
        "CellType",
        drop = TRUE
    ]
    object
}
