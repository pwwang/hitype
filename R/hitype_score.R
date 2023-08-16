#' Calculate ScType scores and assign cell types
#'
#' GNU General Public License v3.0
#'
#' @author
#' Matt Mulvahill, Panwen Wang
#'
#' @importFrom dplyr %>%
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
hitype_score <- function(exprs, gs, scaled = FALSE){
    # Check input
    if (is.data.frame(exprs)) {
        exprs <- as.matrix(exprs)
    }
    if (!is.matrix(exprs)) {
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

    all_markers <- unlist(
        lapply(gs, function(x) unlist(lapply(x, function(y) y$markers)))
    )
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

    # Filter the genes in gs
    gs <- lapply(gs, function(x) {
        lapply(x, function(y) {
            y$markers <- intersect(y$markers, all_markers)
            y$weights <- y$weights[which(y$markers %in% all_markers)]
            y
        })
    })

    lapply(
        gs,
        function(gs_level) hitype_score_level(exprs, gs_level, scaled)
    )
}

#' Calculate ScType scores and assign cell types for one level
#' @param exprs Input scRNA-seq data matrix
#'  (rownames - gene names, colnames - cell names)
#' @param gs_level One level of gene sets prepared by \code{\link{gs_prepare}}
#' @param scaled Whether the input data is scaled or not
#' @return A matrix of cell type scores for each cell
#'  (rownames - cell types, colnames - cell names)
hitype_score_level <- function(exprs, gs_level, scaled = FALSE) {
    # gs_level
    # list(
    #     CD4 = list(
    #       # CD4++++, IL2RA+++, IL2RB++, IL2RG+, IL7R
    #       markers = c("CD4", "IL2RA", "IL2RB", "IL2RG", "IL7R", "CD68"),
    #       weights = c(1, .8, .6, .4, .2, -1)
    #     ),
    #     CD8 = list(
    #       # CD8A++++, CD8B++++, IL2RA+++, IL2RB++, IL2RG+
    #       markers = c("CD8A", "CD8B", "IL2RA", "IL2RB", "IL2RG", "CD68"),
    #       weights = c(1, 1, .8, .6, .4, -1)
    #     )
    # )
    all_markers <- unlist(
        lapply(
            gs_level,
            function(x) x$markers
        )
    )

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

    Z <- if (scaled) exprs else t(scale(t(exprs)))
    Z <- Z[marker_sensitivity$gene_, ]

    # multiple by marker sensitivity
    Z <- Z * marker_sensitivity$score_marker_sensitivity

    gfun <- function(gss_) {
        #       cell1 cell2 cell3
        # Gene1   1     2     3
        # Gene2   4     5     6
        gs_z = Z[
            gs_level[[gss_]]$markers, , drop=FALSE
        ] * gs_level[[gss_]]$weights
        len_z = sqrt(nrow(gs_z))
        colSums(gs_z) / len_z
    }

    es = data.frame(t(
        matrix(
            unlist(lapply(names(gs_level), gfun)),
            ncol = length(gs_level)
        ))
    )

    dimnames(es) = list(names(gs_level), colnames(Z))
    es_max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows

    es_max
}

EMPTY <- "<EMPTY>"  # nolint
UNKNOWN <- "<UNKNOWN>"  # nolint

#' Assign cell types based on ScType scores
#'
#' @param clusters A named vector of original cluster assignments
#'  (names - cell names, values - cluster assignments)
#' @param hitype_scores A list of matrices of cell type scores for each level
#' @param threshold A threshold for low confidence cell type assignment
#'  The cell types are only assigned for cells with scores higher than the
#'  threshold * ncells.
#'  (0 - 1, default 0.25)
#' @param top The number of top cell types to assign for each cluster
#' @return A data from of top cell type assignments
hitype_assign_level <- function(
    clusters,
    hitype_scores,
    threshold,
    top = 10
) {
    x <- do.call(
        "rbind",
        lapply(unique(clusters), function(cl) {
            es_max_cl <- sort(
                rowSums(hitype_scores[ , names(clusters[clusters == cl])]),
                decreasing = TRUE
            )
            head(
                data.frame(
                    Cluster = cl,
                    CellType = names(es_max_cl),
                    Scores = es_max_cl,
                    NCells = sum(clusters == cl)
                ),
                top
            )
        })
    )
    x %>% dplyr::mutate(
        CellType = dplyr::if_else(
            Scores < NCells * threshold,
            UNKNOWN,
            CellType
        )
    )
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
                return (NULL)
            }
        } else {
            if (length(gs) == 1 && gs == EMPTY) {
                return (out)
            }
            if (ct %in% gs) {
                if (ct != UNKNOWN) { out <- paste(c(out, ct), collapse = " ") }
                break
            } else {
                return (NULL)
            }
        }
    }
    return(out)
}

#' Assign cell types based on ScType scores
#'
#' @param clusters A named vector of original cluster assignments
#'  (names - cell names, values - cluster assignments)
#' @param hitype_scores A list of matrices of cell type scores for each level
#' @param gs The gene sets prepared by \code{\link{gs_prepare}}
#'  The `cell_names` is actually used. One could also pass `gs$cell_names`
#'  directly.
#' @param fallback A fallback cell type if no cell type is assigned
#' @param threshold A threshold for low confidence cell type assignment
#'  The cell types are only assigned for cells with scores higher than the
#'  threshold * ncells.
#'  (0 - 1, default 0.25)
#' @return A list of mappings from original cluster assignments to cell types
#'
#' @export
hitype_assign <- function(
    clusters,
    hitype_scores,
    gs = NULL,
    fallback = "Unknown",
    threshold = 0.05
) {
    uclusters <- unique(clusters)
    if (is.data.frame(hitype_scores) || is.matrix(hitype_scores)) {
        # to be compatible with sctype
        cl_results <- hitype_assign_level(
            clusters,
            hitype_scores,
            threshold,
            top = 2
        )

        scores <- cl_results %>%
            dplyr::group_by(Cluster) %>%
            dplyr::slice_max(Scores, n = 1, with_ties = TRUE)

        if (nrow(scores) > length(uclusters)) {
            hitype_scores_count <- scores %>%
                dplyr::count(Cluster) %>%
                dplyr::filter(n > 1)
            warning(
                paste0(
                    paste(
                        "Scores tied in the following clusters:",
                        paste(hitype_scores_count$Cluster, collapse = ", ")
                    ),
                    "\nScores matrix:\n",
                    paste(capture.output(scores), collapse = "\n")
                ),
                immediate. = TRUE
            )
        }

        return(
            scores %>%
                dplyr::slice_max(Scores, n = 1, with_ties = FALSE) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(
                    CellType = dplyr::if_else(
                        CellType == UNKNOWN,
                        fallback,
                        CellType
                    ),
                    CellType = make.unique(CellType)
                )
        )
    }

    if (is.null(gs)) {
        stop("hitype_assign: `gs` is required for multi-level assignment.")
    }
    if ("cell_names" %in% names(gs)) {
        gs <- gs$cell_names
    }

    # Cluster, CellType, Scores, NCells, Level
    cl_results <- do.call(
        rbind,
        lapply(seq_along(hitype_scores), function(i) {
            cl_ret <- hitype_assign_level(
                clusters,
                hitype_scores[[i]],
                threshold
            )
            cl_ret$Level <- i
            cl_ret
        })
    )
    # Work on each Cluster
    cl_results <- lapply(
        split(cl_results, cl_results$Cluster),
        function(cl_ret) {
            cell_types <- lapply(unique(cl_ret$Level), function(l) {
                cl_ret[cl_ret$Level == l, "CellType"]
            })
            # Make sure first level listed in order
            # 1 1 2 2 instead of 1 2 1 2
            cell_types <- rev(cell_types)
            all_types <- do.call(expand.grid, cell_types)
            all_types <- all_types[, rev(colnames(all_types)), drop = FALSE]
            if (nrow(all_types) == 0) {
                return(fallback)
            }
            for (i in seq_len(nrow(all_types))) {
                row <- as.character(unname(unlist(all_types[i, , drop = TRUE])))
                ct <- valid_cell_type(row, gs)
                if (!is.null(ct)) {
                    return(ct)
                }
            }
            return(fallback)
        }
    )
    data.frame(Cluster = names(cl_results), CellType = unlist(cl_results))
}