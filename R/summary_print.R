#' Summarize the hitype_result object
#'
#' @importFrom dplyr %>%
#'
#' @param top The top assigned cell types to return for each cluster
#' @param level_weights The weights for each level of the hierarchy to calculate
#'  the final cell type score
#'  It should be either a numeric vector of length equal to the number of levels
#'  or a single numeric value to be used for all levels
#' @param ... Additional arguments to pass to \code{\link{print}}
#'
#' @return The summary of the hitype_result object
#'
#' @export
summary.hitype_result <- function(hitype_res, top = 1, level_weights = 1) {
    ulevels <- unique(hitype_res$Level)
    if (length(ulevels) == 1) {
        hitype_res %>%
            dplyr::group_by(Cluster) %>%
            dplyr::slice_max(Score, n = top) %>%
            dplyr::ungroup()
    } else {
        if (length(level_weights) == 1) {
            level_weights <- rep(level_weights, length(ulevels))
        }

        # Work on each Cluster
        cl_results <- lapply(
            split(hitype_res, hitype_res$Cluster),
            function(cl_ret) {
                cell_types <- lapply(ulevels, function(l) {
                    cl_ret[cl_ret$Level == l, "CellType"]
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
                        CellType = fallback,
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
                        "Score"
                    ]
                    all_scores <- all_scores + scores / (10^(i - 1))
                }
                all_types$Score <- all_scores
                all_types$Cluster <- cl_ret$Cluster[1]
                all_types <- all_types[order(-all_scores), , drop = FALSE]
                all_types <- all_types %>%
                    dplyr::rowwise() %>%
                    dplyr::mutate(
                        CellType = valid_cell_type(
                            dplyr::c_across(dplyr::starts_with("Level")),
                            attr(hitype_res, "gs")
                        )
                    ) %>%
                    dplyr::ungroup() %>%
                    dplyr::filter(!is.na(CellType)) %>%
                    dplyr::select(Cluster, CellType, Score) %>%
                    dplyr::slice_max(Score, n = top)

                if (nrow(all_types) == 0) {
                    # nocov start
                    return(data.frame(
                        Cluster = cl_ret$Cluster[1],
                        CellType = fallback,
                        Score = NA
                    ))
                    # nocov end
                }
                head(all_types, top)
            }
        )
        do_call(rbind, cl_results) %>%
            dplyr::mutate(CellType = make.unique(CellType)) %>%
            dplyr::arrange(Cluster, dplyr::desc(Score))
    }
}

#' Print the summary of the hitype_result object
#'
#' @param top The top assigned cell types to return for each cluster
#' @param level_weights The weights for each level of the hierarchy to calculate
#'  the final cell type score
#'  It should be either a numeric vector of length equal to the number of levels
#'  or a single numeric value to be used for all levels
#' @param ... Additional arguments to pass to \code{\link{print}}
#'
#' @return The summary of the hitype_result object
#'
#' @export
print.hitype_result <- function(hitype_res, top = 1, level_weights = 1, ...) {
    # nocov start
    x <- summary(hitype_res, top = top, level_weights = level_weights)
    print(x, ...)
    invisible(x)
    # nocov end
}

#' Get the valid cell type for a given cell type
#'
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
