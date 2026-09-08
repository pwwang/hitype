#' Summarize the hitype_result object
#'
#' @importFrom dplyr %>%
#' @importFrom utils head
#'
#' @param object A hitype_result object
#' @param top The top assigned cell types to return for each cluster
#' @param level_weights The weights for each level of the hierarchy to calculate
#'  the final cell type score
#'  It should be either a numeric vector of length equal to the number of levels
#'  or a single numeric value to be used for all levels
#'  It can also be a function that takes the levels as input and returns a
#'  numeric vectors as the weights.
#' @param make_unique Whether to make the cell type names unique
#' @param ... Additional arguments passed to the specific method.
#'
#' @return The summary of the hitype_result object
#'
#' @export
summary.hitype_result <- function(
    object,
    top = 1,
    level_weights = function(l) 1 / (2 ^ (l - 1)),
    make_unique = FALSE,
    ...
) {
    ulevels <- unique(object$Level)
    if (length(ulevels) == 1) {
        out <- object %>%
            dplyr::group_by(Cluster) %>%
            dplyr::slice_max(Score, n = top, with_ties = FALSE) %>%
            dplyr::ungroup()
    } else {
        if (is.function(level_weights)) {
            level_weights <- level_weights(ulevels)
        }
        if (length(level_weights) == 1) {
            level_weights <- rep(
                level_weights / length(ulevels),
                length(ulevels)
            )
        }

        # Work on each Cluster
        cl_results <- lapply(
            split(object, object$Cluster),
            function(cl_ret) {
                cell_types <- lapply(ulevels, function(l) {
                    cl_ret[cl_ret$Level == l, "CellType", drop = TRUE]
                })
                # Make sure first level listed in order
                # 1 1 2 2 instead of 1 2 1 2
                cell_types <- rev(cell_types)
                cell_types$stringsAsFactors <- FALSE
                all_types <- do_call(expand.grid, cell_types)
                all_types <- all_types[, rev(colnames(all_types)), drop = FALSE]
                if (nrow(all_types) == 0) {
                    # nocov start
                    return(data.frame(
                        Cluster = cl_ret$Cluster[1],
                        CellType = attr(object, "fallback"),
                        Score = NA
                    ))
                    # nocov end
                }
                colnames(all_types) <- paste0("Level", ulevels)
                # Add Scores to all_types
                all_scores <- rep(0, nrow(all_types))
                for (i in seq_along(ulevels)) {
                    scores <- cl_ret[cl_ret$Level == ulevels[i], , drop = FALSE]
                    scores <- scores[
                        match(all_types[, i], scores$CellType),
                        "Score",
                        drop = TRUE
                    ]
                    all_scores <- all_scores + scores * level_weights[i]
                }
                all_types$Score <- all_scores
                all_types$Cluster <- cl_ret$Cluster[1]
                all_types <- all_types[order(-all_scores), , drop = FALSE]
                all_types <- all_types %>%
                    dplyr::rowwise() %>%
                    dplyr::mutate(
                        CellType = valid_cell_type(
                            dplyr::c_across(dplyr::starts_with("Level")),
                            attr(object, "gs")
                        )
                    ) %>%
                    dplyr::ungroup() %>%
                    dplyr::filter(!is.na(CellType)) %>%
                    dplyr::select(Cluster, CellType, Score) %>%
                    dplyr::distinct(Cluster, CellType, .keep_all = TRUE) %>%
                    dplyr::slice_max(Score, n = top, with_ties = FALSE)

                if (nrow(all_types) == 0) {
                    # nocov start
                    return(data.frame(
                        Cluster = cl_ret$Cluster[1],
                        CellType = attr(object, "fallback"),
                        Score = NA
                    ))
                    # nocov end
                }
                head(all_types, top)
            }
        )
        out <- do_call(rbind, cl_results) %>%
            # dplyr::mutate(CellType = make.unique(CellType)) %>%
            dplyr::arrange(Cluster, dplyr::desc(Score))
    }
    if (make_unique) {
        out$CellType <- make.unique(out$CellType)
    }
    out
}

#' Print the summary of the hitype_result object
#'
#' @param x A hitype_result object
#' @param top The top assigned cell types to return for each cluster
#' @param level_weights The weights for each level of the hierarchy to calculate
#'  the final cell type score
#'  It should be either a numeric vector of length equal to the number of levels
#'  or a single numeric value to be used for all levels
#'  It can also be a function that takes the levels as input and returns a
#'  numeric vectors as the weights.
#' @param make_unique Whether to make the cell type names unique
#' @param ... Additional arguments to pass to \code{\link{print}}
#'
#' @return The summary of the hitype_result object
#'
#' @export
print.hitype_result <- function(
    x,
    top = 1,
    level_weights = function(l) 1 / (2 ^ (l - 1)),
    make_unique = FALSE,
    ...
) {
    s <- summary(x, top = top, level_weights = level_weights,
        make_unique = make_unique)
    print(s, ...)
    invisible(s)
}

#' Get the valid cell type for a given cell type
#'
#' @keywords internal
#' @param cell_type A vector cell types, with length the same as the number of
#'  cell type levels
#' @param gs A list of cell names for each cell type level
#' @return A string of valid cell types or NULL if the cell type is invalid
valid_cell_type <- function(cell_type, gs) {
    out <- NULL
    for (ct in cell_type) {
        if (is.list(gs)) {
            if (ct %in% names(gs)) {
                if (ct != UNKNOWN) { out <- paste(c(out, ct), collapse = " ") }
                gs <- gs[[ct]]
            } else {
                return(NA_character_)
            }
        } else {
            if (length(gs) == 1 && gs == EMPTY) {
                return(out)
            }
            if (ct %in% gs) {
                if (ct != UNKNOWN) { out <- paste(c(out, ct), collapse = " ") }
                break
            } else {
                return(NA_character_)
            }
        }
    }
    if (is.null(out)) NA_character_ else out
}
