# GNU General Public License v3.0

#' Calculate cell type scores
#'
#' @importFrom dplyr %>%
#'
#' @author Matt Mulvahill, Panwen Wang
#'
#' @param exprs Input scRNA-seq data matrix
#'  (rownames - gene names, colnames - cell names)
#' @param gs The gene sets prepared by \code{\link{gs_prepare}}
#'  The `gene_sets` is used. One could also pass `gs$gene_sets` directly.
#' @param scaled Whether the input data is scaled or not
#'
#' @return A list of matrices of cell type scores for each level
#'  (rownames - cell types, colnames - cell names)
#' @export
hitype_score <- function(exprs, gs, scaled = FALSE) {
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
        function(gs_level) hitype_score_level(z, gs_level)
    )
}

#' Calculate ScType scores and assign cell types for one level
#'
#' @keywords internal
#'
#' @param z Z-scaled expression matrix
#'  (rownames - gene names, colnames - cell names)
#' @param gs_level One level of gene sets prepared by \code{\link{gs_prepare}}
#'
#' @return A matrix of cell type scores for each cell
#'  (rownames - cell types, colnames - cell names)
hitype_score_level <- function(z, gs_level) {
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
    # if (!has_plus_minus) {
        # multiple by marker sensitivity
        z <- z * marker_sensitivity$score_marker_sensitivity
    # }

    gfun <- function(gss_) {
        #       cell1 cell2 cell3
        # Gene1   1     2     3
        # Gene2   4     5     6
        gs_z <- z[
            gs_level[[gss_]]$markers, , drop = FALSE
        ] * gs_level[[gss_]]$weights
        colSums(gs_z) / sqrt(nrow(gs_z))
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
#' @param scores A list of matrices of cell type scores for each level
#' @param threshold A threshold for low confidence cell type assignment
#'  The cell types are only assigned for cells with scores higher than the
#'  threshold * ncells.
#'  (0 - 1, default 0.05)
#' @param top The number of top cell types to assign for each cluster
#' @return A data from of top cell type assignments
hitype_assign_level <- function(
    clusters,
    scores,
    threshold,
    top = 10
) {
    scores <- scales::rescale(as.matrix(scores), to = c(0, 1))
    x <- do_call(
        "rbind",
        lapply(unique(clusters), function(cl) {
            ncells <- sum(clusters == cl)
            es_max_cl <- sort(
                rowSums(
                    scores[, names(clusters[clusters == cl]), drop = FALSE]
                ) / ncells,
                decreasing = TRUE
            )
            head(
                data.frame(
                    Cluster = cl,
                    CellType = names(es_max_cl),
                    Score = es_max_cl
                ),
                top
            )
        })
    )
    x %>% dplyr::mutate(
        CellType = dplyr::if_else(Score < threshold, UNKNOWN, CellType)
    )
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
#' @param threshold A threshold for low confidence cell type assignment
#'  The cell types are only assigned for cells with scores higher than the
#'  threshold * ncells.
#'  (0 - 1, default 0.05)
#' @param top The number of top cell types to assign for each cluster in the
#'  result.
#' @return A dataframe with columns: `Level`, `Cluster`, `CellType` and `Score`.
#'  For each level and cluster, the top cell types are returned.
#'  You can use \code{\link{summary.hitype_result}} to print the combination of
#'  cell types for each cluster.
#'
#' @export
hitype_assign <- function(
    clusters,
    scores,
    gs = NULL,
    fallback = "Unknown",
    threshold = 0.05,
    top = 10
) {
    if (is.data.frame(scores) || is.matrix(scores)) {
        # to be compatible with sctype
        # Cluster, CellType, Score
        out <- hitype_assign_level(clusters, scores, threshold, top = top) %>%
            dplyr::mutate(
                Level = 1,
                CellType = dplyr::if_else(
                    CellType == UNKNOWN,
                    fallback,
                    CellType
                )
            ) %>%
            dplyr::select(Level, Cluster, CellType, Score)
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
                    clusters, scores[[i]], threshold, top = top
                )
                cl_ret$Level <- i
                cl_ret
            })
        )
    }

    cluster_order <- suppressWarnings(as.numeric(out$Cluster))
    out <- out[
        order(as.numeric(out$Level), cluster_order, -out$Score),
        c("Level", "Cluster", "CellType", "Score"),
        drop = FALSE
    ]
    rownames(out) <- NULL
    class(out) <- c("hitype_result", class(out))
    attr(out, "gs") <- gs
    attr(out, "fallback") <- fallback
    return(out)
}