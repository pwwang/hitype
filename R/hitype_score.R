# GNU General Public License v3.0

#' Calculate cell type scores
#'
#' @importFrom dplyr %>%
#' @importFrom stats na.omit
#' @importFrom utils head
#'
#' @author Matt Mulvahill, Panwen Wang
#'
#' @param exprs Input scRNA-seq data matrix
#'  (rownames - gene names, colnames - cell names)
#' @param gs The gene sets prepared by \code{\link{gs_prepare}}
#'  The `gene_sets` is used. One could also pass `gs$gene_sets` directly.
#' @param scaled Whether the input data is scaled or not
#' @param norm Normalization method for scoring. \code{"sqrt"} (default)
#'  divides by \code{sqrt(n)}; \code{"weight"} divides by
#'  \code{sum(abs(weights))}.
#' @param use_sensitivity Whether to multiply expression z-scores by
#'  marker sensitivity scores. Set to \code{FALSE} when using learned
#'  weights to avoid double-penalizing shared markers. Default is \code{TRUE}.
#'
#' @return A list of matrices of cell type scores for each level
#'  (rownames - cell types, colnames - cell names)
#' @export
hitype_score <- function(exprs, gs, scaled = FALSE, norm = "sqrt",
                        use_sensitivity = TRUE) {
    # Check input
    if (is.data.frame(exprs)) {
        exprs <- as.matrix(exprs)
    }
    if (
        !is.matrix(exprs) &&
        !any(class(exprs) %in% c("dgCMatrix", "dgTMatrix"))
    ) {
        stop("Input scRNA-seq data must be a matrix or data.frame")
    }
    if (sum(dim(exprs) == 0)) {
        stop("Input scRNA-seq data is empty")
    }
    if (!is.list(gs)) {
        stop("Input gene sets must be a list")
    }
    if ("gene_sets" %in% names(gs)) {
        gs <- gs$gene_sets
    }
    if (length(gs) == 0) {
        stop("Input gene sets is empty")
    }
    if (!is.logical(scaled)) {
        stop("Input scaled must be a logical value")
    }

    # gs
    # list(
    #   # Level 1
    #   list(
    #       CD4 = list(
    #         # CD4++++, IL2RA+++, IL2RB++, IL2RG+, IL7R
    #         markers = c("CD4", "IL2RA", "IL2RB", "IL2RG", "IL7R", "CD68"),
    #         weights = c(1, .8, .6, .4, .2, -1)
    #       ),
    #       CD8 = list(
    #         # CD8A++++, CD8B++++, IL2RA+++, IL2RB++, IL2RG+
    #         markers = c("CD8A", "CD8B", "IL2RA", "IL2RB", "IL2RG", "CD68"),
    #         weights = c(1, 1, .8, .6, .4, -1)
    #       )
    #   ),
    #   # Level 2
    #   list(
    #       Naive = list(
    #         markers = c("CCR7", "SELL", "CD27"),
    #         weights = c(1, .8, .6)
    #       ),
    #       Memory = list(
    #         markers = c("CD44", "CD69", "CD45RA", "CD45RO"),
    #         weights = c(1, .8, .6, .4)
    #       )
    #   )
    # )

    all_markers <- na.omit(unlist(
        lapply(gs, function(x) unlist(lapply(x, function(y) y$markers)))
    ))
    # Check if all markers are in the input data
    non_exist_markers <- setdiff(all_markers, rownames(exprs))
    if (length(non_exist_markers) > 0) {
        w <- paste0(
            "There are markers not in the input data. Showing:\n",
            "    Orignal Marker -> Suggested Correction -> Suggested in Input\n"
        )
        for (ne_marker in non_exist_markers) {
            x <- suppressMessages(suppressWarnings(
                HGNChelper::checkGeneSymbols(ne_marker)$Suggested.Symbol
            ))
            suggested_in_input <- if (is.na(x)) {
                NA
            } else {
                x %in% rownames(exprs)
            }
            w <- paste0(
                w,
                paste0(
                    "    ",
                    ne_marker, " -> ", x, " -> ", suggested_in_input, "\n"
                )
            )
        }
        warning(w, immediate. = TRUE)
    }
    all_markers <- intersect(all_markers, rownames(exprs))
    if (length(all_markers) == 0) {
        stop("No markers are in the input data. Are they in the same format?")
    }

    exprs <- exprs[all_markers, , drop = FALSE]
    # Filter the genes in gs
    gs <- lapply(gs, function(x) {
        lapply(x, function(y) {
            y$markers <- intersect(y$markers, all_markers)
            y$weights <- y$weights[which(y$markers %in% all_markers)]
            y
        })
    })

    if (scaled) { z <- exprs }  # nocov
    else if (any(class(exprs) %in% c("dgCMatrix", "dgTMatrix"))) {
        z <- t(scale(Matrix::t(exprs)))  # nocov
    } else {
        z <- t(scale(t(exprs)))
    }

    lapply(
        gs,
        function(gs_level) hitype_score_level(z, gs_level, norm, use_sensitivity)
    )
}

#' Calculate ScType scores and assign cell types for one level
#'
#' @keywords internal
#'
#' @param z Z-scaled expression matrix
#'  (rownames - gene names, colnames - cell names)
#' @param gs_level One level of gene sets prepared by \code{\link{gs_prepare}}
#' @param norm Normalization method for scoring. \code{"sqrt"} (default)
#'  divides by \code{sqrt(n)}; \code{"weight"} divides by
#'  \code{sum(abs(weights))}.
#' @param use_sensitivity Whether to multiply expression z-scores by
#'  marker sensitivity scores. Set to \code{FALSE} when using learned
#'  weights to avoid double-penalizing shared markers. Default is \code{TRUE}.
#'
#' @return A matrix of cell type scores for each cell
#'  (rownames - cell types, colnames - cell names)
hitype_score_level <- function(z, gs_level, norm = "sqrt",
                             use_sensitivity = TRUE) {
    # gs_level
    # list(
    #     CD4 = list(
    #       # CD4++++, IL2RA+++, IL2RB++, IL2RG+, IL7R
    #       markers = c("CD4", "IL2RA", "IL2RB", "IL2RG", "IL7R", "CD68"),
    #       weights = c(5, 4, 3, 2, 1)
    #     ),
    #     CD8 = list(
    #       # CD8A++++, CD8B++++, IL2RA+++, IL2RB++, IL2RG+
    #       markers = c("CD8A", "CD8B", "IL2RA", "IL2RB", "IL2RG", "CD68"),
    #       weights = c(5, 5, 4, 3, 2)
    #     )
    # )
    all_markers <- unlist(lapply(gs_level, function(x) x$markers))

    # Marker stat
    # PTPRC    ITGAM    IL2RA    CD247
    # 18       17       16       13
    # ISG20     SELL    IL3RA     CD24
    # 10       10        9        8
    marker_stat <- sort(table(all_markers), decreasing = T)
    # Marker sensitivity
    #     score_marker_sensitivity    gene_
    # 1                        -16    PTPRC
    # 2                        -15    ITGAM
    # 3                        -14    IL2RA
    # 4                        -11    CD247
    # 5                        -11     CD3D
    # 6                        -11     CD3E
    # 7                        -11     CD3G
    # 8                         -9     CD27
    marker_sensitivity <- data.frame(
        score_marker_sensitivity = scales::rescale(
            as.numeric(marker_stat),
            to = c(0, 1),
            from = c(length(gs_level), 1)
        ),
        gene_ = names(marker_stat),
        stringsAsFactors = FALSE
    )

    z <- z[marker_sensitivity$gene_, ]

    # has_plus_minus <- FALSE
    # for (gs_ in gs_level) {
    #     if (any(grepl("[+-]$", gs_$markers))) {
    #         has_plus_minus <- TRUE
    #         break
    #     }
    # }
    # Multiple by marker sensitivity (can be disabled for learned weights)
    if (use_sensitivity) {
        z <- z * marker_sensitivity$score_marker_sensitivity
    }

    gfun <- function(gss_) {
        #       cell1 cell2 cell3
        # Gene1   1     2     3
        # Gene2   4     5     6
        gs_z <- z[
            gs_level[[gss_]]$markers, , drop = FALSE
        ] * gs_level[[gss_]]$weights
        if (norm == "weight") {
            w_sum <- sum(abs(gs_level[[gss_]]$weights))
            if (w_sum == 0) w_sum <- 1
            colSums(gs_z) / w_sum
        } else {
            colSums(gs_z) / sqrt(nrow(gs_z))
        }
    }

    es = data.frame(t(
        matrix(
            unlist(lapply(names(gs_level), gfun)),
            ncol = length(gs_level)
        ))
    )

    dimnames(es) <- list(names(gs_level), colnames(z))
    es_max <- es[!apply(is.na(es) | es == "", 1, all), ] # remove na rows

    es_max
}

#' Assign cell types based on ScType scores
#'
#' @keywords internal
#'
#' @param clusters A named vector of original cluster assignments
#'  (names - cell names, values - cluster assignments)
#' @param scores A matrix of cell type scores (cell_types x cells)
#' @param threshold Confidence threshold as top1/top2 score ratio.
#'  When `NULL` (default), no filtering is done.
#'  When a number, the rank-1 cell type of a cluster is marked as
#'  `<UNKNOWN>` when its score is less than `threshold` times the
#'  second-best score.
#' @param top The number of top cell types to assign for each cluster
#' @param mode `"cluster"` (default) aggregates scores per cluster
#'  then assigns. `"cell"` assigns each cell individually then
#'  reports majority vote per cluster.
#' @return A data from of top cell type assignments with columns:
#'  Cluster, CellType, Score, Margin
hitype_assign_level <- function(
    clusters,
    scores,
    threshold = NULL,
    top = 10,
    mode = c("cluster", "cell")
) {
    mode <- match.arg(mode)
    scores <- as.matrix(scores)

    if (mode == "cell") {
        # Per-cell assignment: top-scoring cell type per cell
        cell_types <- rownames(scores)
        cell_assignments <- do_call("rbind", lapply(
            seq_len(ncol(scores)),
            function(j) {
                col <- scores[, j]
                ord <- order(col, decreasing = TRUE)
                data.frame(
                    best = cell_types[ord[1]],
                    score1 = col[ord[1]],
                    best2 = cell_types[ord[2]],
                    score2 = col[ord[2]],
                    stringsAsFactors = FALSE
                )
            }
        ))
        cell_assignments$Cluster <- clusters[colnames(scores)]
        cell_assignments$Margin <- cell_assignments$score1 -
            cell_assignments$score2

        # Majority vote per cluster
        x <- do_call("rbind", lapply(unique(clusters), function(cl) {
            cl_cells <- cell_assignments[
                cell_assignments$Cluster == cl, , drop = FALSE]
            ncells <- nrow(cl_cells)
            type_counts <- sort(table(cl_cells$best), decreasing = TRUE)
            total_margin <- tapply(cl_cells$Margin,
                cl_cells$best, sum)
            data.frame(
                Cluster = cl,
                CellType = names(type_counts),
                Score = as.numeric(type_counts) / ncells,
                Margin = as.numeric(total_margin[names(type_counts)]),
                stringsAsFactors = FALSE
            )
        }))
    } else {
        # Cluster mode: mean score per cell type within each cluster
        x <- do_call(
            "rbind",
            lapply(unique(clusters), function(cl) {
                ncells <- sum(clusters == cl)
                cl_scores <- rowSums(
                    scores[, names(clusters[clusters == cl]),
                           drop = FALSE]
                ) / ncells
                ord <- order(cl_scores, decreasing = TRUE)
                head(
                    data.frame(
                        Cluster = cl,
                        CellType = names(cl_scores)[ord],
                        Score = cl_scores[ord],
                        Margin = c(
                            diff(-cl_scores[ord]),
                            NA
                        )[seq_along(ord)],
                        stringsAsFactors = FALSE
                    ),
                    top
                )
            })
        )
    }

    # Apply threshold as top1/top2 score ratio
    if (!is.null(threshold) && is.numeric(threshold)) {
        x <- x %>%
            dplyr::group_by(Cluster) %>%
            dplyr::mutate(
                CellType = dplyr::if_else(
                    dplyr::row_number() == 1 &
                        Score < threshold *
                        dplyr::lead(Score, default = Inf),
                    UNKNOWN, CellType
                )
            ) %>%
            dplyr::ungroup()
    }

    x
}

#' Generate scores for cell types for each level
#'
#' @param clusters A named vector of original cluster assignments
#'  (names - cell names, values - cluster assignments)
#' @param scores A list of matrices of cell type scores for each level
#' @param gs The gene sets prepared by \code{\link{gs_prepare}}
#'  The `cell_names` is actually used. One could also pass `gs$cell_names`
#'  directly.
#' @param fallback A fallback cell type if no cell type is assigned
#' @param threshold Confidence threshold as top1/top2 score ratio.
#'  `NULL` (default) means no confidence filtering.
#'  A number marks the rank-1 cell type of a cluster as `<UNKNOWN>`
#'  when its score is less than `threshold` times the second-best score.
#' @param top The number of top cell types to assign for each cluster in the
#'  result.
#' @param mode `"cluster"` (default) aggregates scores per cluster
#'  then assigns the top cell type. `"cell"` assigns each cell
#'  individually then reports majority vote per cluster.
#' @return A dataframe with columns: `Level`, `Cluster`, `CellType`, `Score`,
#'  and `Margin`. For each level and cluster, the top cell types are returned.
#'  You can use \code{\link{summary.hitype_result}} to print the combination of
#'  cell types for each cluster.
#'
#' @export
hitype_assign <- function(
    clusters,
    scores,
    gs = NULL,
    fallback = "Unknown",
    threshold = NULL,
    top = 10,
    mode = c("cluster", "cell")
) {
    mode <- match.arg(mode)
    if (is.data.frame(scores) || is.matrix(scores)) {
        # to be compatible with sctype
        # Cluster, CellType, Score
        out <- hitype_assign_level(
            clusters, scores, threshold, top = top, mode = mode
        ) %>%
            dplyr::mutate(
                Level = 1,
                CellType = dplyr::if_else(
                    CellType == UNKNOWN,
                    fallback,
                    CellType
                )
            ) %>%
            dplyr::select(Level, Cluster, CellType, Score, Margin)
    } else {
        if (is.null(gs)) {
            stop("hitype_assign: `gs` is required for multi-level assignment.")
        }
        if ("cell_names" %in% names(gs)) {
            gs <- gs$cell_names
        }

        # Cluster, CellType, Score, Level
        out <- do_call(
            rbind,
            lapply(seq_along(scores), function(i) {
                cl_ret <- hitype_assign_level(
                    clusters, scores[[i]], threshold, top = top, mode = mode
                )
                cl_ret$Level <- i
                cl_ret
            })
        )
    }

    cluster_order <- suppressWarnings(as.numeric(out$Cluster))
    out <- out[
        order(as.numeric(out$Level), cluster_order, -out$Score),
        c("Level", "Cluster", "CellType", "Score", "Margin"),
        drop = FALSE
    ]
    rownames(out) <- NULL
    class(out) <- c("hitype_result", class(out))
    attr(out, "gs") <- gs
    attr(out, "fallback") <- fallback
    return(out)
}