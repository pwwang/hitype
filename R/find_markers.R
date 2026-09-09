#' Find marker genes for cell types
#'
#' @description
#' Discovers marker genes for each cell type (cluster) in an expression
#' matrix or Seurat object, and returns them in the hitype marker database
#' format so the result can be passed directly to [gs_prepare()]. Three
#' backends are available: the dependency-light fold-change method
#' (default), Seurat's `FindAllMarkers`, and presto's Wilcoxon rank-sum
#' test. For `method = "fc"`, the input expression matrix is expected to
#' be log-normalized.
#'
#' @param exprs A genes x cells expression matrix (plain matrix or
#'  dgCMatrix/dgTMatrix) or a Seurat object. For `method = "fc"`, the
#'  input is expected to be log-normalized.
#' @param clusters A named vector of cluster ids (names = cells) or a
#'  factor. If `exprs` is a Seurat object, defaults to
#'  `Seurat::Idents(exprs)`.
#' @param method The marker-finding method. One of:
#'  \describe{
#'    \item{"fc"}{Fold-change based (default). No extra dependencies.
#'      Ranks genes by `log2fc * (pct_in - pct_out)`.}
#'    \item{"seurat"}{Seurat's `FindAllMarkers`. Requires a Seurat
#'      object.}
#'    \item{"presto"}{presto's Wilcoxon rank-sum test (`wilcoxauc`).
#'      Requires the presto package.}
#'  }
#' @param top Number of markers to return per cell type: a single number
#'  applied to both positive and negative markers, or a length-2 vector
#'  `c(n_positive, n_negative)` budgeting each direction separately
#'  (default `c(10, 10)`). The negative budget is only used when
#'  `pos_only = FALSE`.
#' @param min_log2fc Minimum log2 fold change for a gene to be kept as a
#'  marker (when `pos_only = TRUE`).
#' @param min_pct Minimum fraction of cells in the cell type expressing
#'  the gene.
#' @param against Restrict the reference group a cell type is compared
#'  against (by default one-vs-rest, i.e. all other cells). One of:
#'  \describe{
#'    \item{`NULL`}{One-vs-rest (default).}
#'    \item{a character vector of cell types}{Only the cells of those
#'      types form the reference group: `pct_out`, `log2fc` and the score
#'      are computed against them instead of all the other cells, which
#'      recovers sibling-specific markers shared with the rest of the
#'      cells. A type listed in `against` cannot use its own cells as the
#'      reference, so they are excluded from its reference group; if no
#'      other listed type remains, no markers are returned for it.}
#'    \item{`"nearest"`}{For each cell type, the single most similar other
#'      type is used as its reference group. Similarity is the Pearson
#'      correlation between the mean expression profiles of the two types
#'      (computed once per call over the input matrix), and the resolved
#'      type is then used exactly like an explicit `against` entry.}
#'  }
#'  An error is raised if `against` names a cell type not present in the
#'  data (available types are listed) or if it names the only cell type
#'  in the data (a type cannot be compared against itself). Negative
#'  markers (`pos_only = FALSE`) are extracted as before as genes low in
#'  the type vs the rest of the cells, but when `against` is set they
#'  must ALSO satisfy `log2fc <= -min_log2fc` against the reference
#'  group, so that sibling-specific negative markers are found.
#'  For `method = "seurat"`, `against` is implemented with per-type
#'  two-group `Seurat::FindMarkers(ident.1 = type, ident.2 = against)`
#'  calls (`"nearest"` is resolved to the per-type reference types
#'  first); `method = "presto"` does not support `against` and raises an
#'  error.
#' @param max_pct_out Drop positive-marker candidates expressed in more
#'  than this fraction of ALL other cells (the one-vs-rest `pct_out`,
#'  computed over all cells not of the type even when `against` is set),
#'  guarding against pan-lineage genes up-regulated in most other cell
#'  types. Negative markers are not subject to this guard. A single
#'  aggregated warning reports the number of dropped candidates per cell
#'  type. Default `0.75`.
#' @param pos_only Only keep genes that are higher in the cell type than
#'  in the rest of the cells (positive markers). If `FALSE`, all genes
#'  passing `min_pct` are considered, ranked by score.
#' @param level The hierarchy level to write into the output data frame.
#' @param format The format of the output data frame. One of
#'   `"universal"` (default) or `"db"` (the hitype/ScType wide format).
#'
#' @return A data frame in the universal marker format (default) with
#'  columns `cell_type`, `gene`, `direction` and `level`, or in the hitype
#'  db format (`format = "db"`) with columns `cellName`,
#'  `geneSymbolmore1`, `geneSymbolmore2` and `level`. Both are directly
#'  consumable by [gs_prepare()].
#'
#' @importFrom stats setNames
#'
#' @examples
#' set.seed(1)
#' ngenes <- 40
#' ncells <- 60
#' exprs <- matrix(runif(ngenes * ncells, 0, 0.2), ngenes, ncells)
#' rownames(exprs) <- paste0("gene", seq_len(ngenes))
#' colnames(exprs) <- paste0("cell", seq_len(ncells))
#' clusters <- setNames(
#'     rep(c("Tcell", "Bcell", "Monocyte"), each = 20),
#'     colnames(exprs)
#' )
#' markers <- find_markers(exprs, clusters, method = "fc", top = 5)
#' head(markers)
#' @export
find_markers <- function(
    exprs,
    clusters = NULL,
    method = c("fc", "seurat", "presto"),
    top = c(10, 10),
    min_log2fc = 0.25,
    min_pct = 0.1,
    against = NULL,
    max_pct_out = 0.75,
    pos_only = TRUE,
    level = 1,
    format = c("universal", "db")
) {
    method <- match.arg(method)
    format <- match.arg(format)
    if (!is.numeric(top) || !length(top) %in% c(1, 2) || anyNA(top) ||
        any(top < 1) || any(top != as.integer(top))) {
        stop(
            "`top` must be a positive whole number (a scalar applies to ",
            "both positive and negative markers) or a length-2 vector ",
            "c(n_positive, n_negative) of positive whole numbers"
        )
    }
    top <- as.integer(rep(top, length.out = 2))
    if (!is.numeric(min_log2fc) || length(min_log2fc) != 1 ||
        is.na(min_log2fc) || min_log2fc < 0) {
        stop("`min_log2fc` must be a single number >= 0")
    }
    if (!is.numeric(min_pct) || length(min_pct) != 1 ||
        is.na(min_pct) || min_pct < 0 || min_pct > 1) {
        stop("`min_pct` must be a single number in [0, 1]")
    }
    if (!is.numeric(max_pct_out) || length(max_pct_out) != 1 ||
        is.na(max_pct_out) || max_pct_out <= 0 || max_pct_out > 1) {
        stop("`max_pct_out` must be a single number in (0, 1]")
    }
    if (!is.null(against)) {
        if (!is.character(against) || length(against) == 0 ||
            anyNA(against)) {
            stop(
                "`against` must be NULL, \"nearest\", or a character ",
                "vector of cell types"
            )
        }
        if (method == "presto") {
            stop(
                "method = \"presto\" does not support `against`; use ",
                "method = \"fc\" (default) or method = \"seurat\""
            )
        }
    }

    is_seurat <- inherits(exprs, "Seurat")
    if (is_seurat) {
        if (is.null(clusters)) {
            clusters <- Seurat::Idents(exprs)
        }
        if (method != "seurat") {
            exprs <- Seurat::GetAssayData(exprs, layer = "data")
        }
    } else {
        if (inherits(exprs, "dgTMatrix")) {
            exprs <- Matrix::Matrix(exprs, sparse = TRUE)
        }
        if (!(is.matrix(exprs) || inherits(exprs, "dgCMatrix"))) {
            stop(
                "`exprs` must be a genes x cells expression matrix ",
                "(plain or dgCMatrix/dgTMatrix) or a Seurat object"
            )
        }
    }

    if (is.null(rownames(exprs))) {
        stop("`exprs` must have rownames (gene symbols)")
    }
    if (nrow(exprs) == 0 || ncol(exprs) == 0) {
        stop("`exprs` must be a non-empty matrix of genes x cells")
    }
    if (is.null(colnames(exprs))) {
        colnames(exprs) <- paste0("cell", seq_len(ncol(exprs)))
    }

    if (is.null(clusters)) {
        stop("Please provide `clusters`: a named vector of cluster ids")
    }
    factor_levels <- NULL
    if (is.factor(clusters)) {
        factor_levels <- levels(clusters)
        cl_names <- names(clusters)
        clusters <- as.character(clusters)
        names(clusters) <- cl_names
    }
    if (!is.null(names(clusters))) {
        if (!setequal(names(clusters), colnames(exprs))) {
            stop("Names of `clusters` must match colnames of `exprs`")
        }
        clusters <- clusters[colnames(exprs)]
    } else {
        if (length(clusters) != ncol(exprs)) {
            stop(
                "`clusters` must have the same length as the number ",
                "of cells (columns) in `exprs`"
            )
        }
        names(clusters) <- colnames(exprs)
    }
    clusters <- as.character(clusters)
    names(clusters) <- colnames(exprs)

    ntab <- table(clusters)
    small <- names(ntab)[ntab < 3]
    if (length(small) > 0) {
        warning(
            "Dropping cell types with fewer than 3 cells: ",
            paste(small, collapse = ", ")
        )
        keep <- !clusters %in% small
        clusters <- clusters[keep]
        exprs <- exprs[, names(clusters)]
    }

    cts <- if (!is.null(factor_levels)) {
        factor_levels[factor_levels %in% unique(clusters)]
    } else {
        unique(clusters)
    }
    if (length(cts) == 0) {
        stop("No cell types left after filtering")
    }
    clusters <- factor(clusters, levels = cts)

    # Reference groups per cell type for `against`. `refs = NULL` keeps
    # the default one-vs-rest comparison; otherwise refs[[ct]] holds the
    # types whose cells form ct's reference group (empty if ct is listed
    # in `against` but no other listed type is left to compare against).
    refs <- NULL
    if (!is.null(against)) {
        if (identical(against, "nearest")) {
            if (length(cts) < 2) {
                stop(
                    "`against = \"nearest\"` requires at least two cell ",
                    "types in the data"
                )
            }
            mat <- if (inherits(exprs, "Seurat")) {
                Seurat::GetAssayData(exprs, layer = "data")
            } else {
                exprs
            }
            # per-type mean expression profiles; correlations between
            # them are computed once per call (cached)
            profs <- setNames(
                lapply(cts, function(ct) {
                    Matrix::rowMeans(mat[, clusters == ct, drop = FALSE])
                }),
                cts
            )
            refs <- setNames(
                lapply(cts, function(ct) {
                    others <- setdiff(cts, ct)
                    cors <- vapply(
                        others,
                        function(other) {
                            stats::cor(profs[[ct]], profs[[other]])
                        },
                        numeric(1)
                    )
                    others[which.max(replace(cors, is.na(cors), -Inf))]
                }),
                cts
            )
        } else {
            missing <- setdiff(against, cts)
            if (length(missing) > 0) {
                stop(
                    "`against` type(s) not present in the data: ",
                    paste(missing, collapse = ", "),
                    ". Available cell types: ",
                    paste(cts, collapse = ", ")
                )
            }
            refs <- setNames(
                lapply(cts, function(ct) setdiff(against, ct)),
                cts
            )
            if (all(lengths(refs) == 0)) {
                stop(
                    "`against` names the type itself: '", cts[1], "' is ",
                    "the only cell type in the data, so it has no other ",
                    "cells to be compared against"
                )
            }
        }
    }

    res <- switch(
        method,
        fc = find_markers_fc(
            exprs, clusters, top, min_log2fc, min_pct, pos_only,
            refs, max_pct_out
        ),
        seurat = find_markers_seurat(
            exprs, clusters, top, min_log2fc, min_pct, pos_only,
            refs, max_pct_out
        ),
        presto = find_markers_presto(
            exprs, clusters, top, min_log2fc, min_pct, pos_only,
            max_pct_out
        )
    )

    db <- data.frame(
        cellName = cts,
        geneSymbolmore1 = unname(res$markers1[cts]),
        geneSymbolmore2 = unname(res$markers2[cts]),
        level = rep(as.integer(level), length(cts)),
        stringsAsFactors = FALSE
    )
    if (format == "universal") {
        parts <- lapply(seq_len(nrow(db)), function(i) {
            p1 <- explode(db$geneSymbolmore1[i]); p1 <- p1[p1 != ""]
            p2 <- explode(db$geneSymbolmore2[i]); p2 <- p2[p2 != ""]
            if (length(p1) == 0 && length(p2) == 0) {
                return(NULL)  # no markers for this cell type
            }
            data.frame(
                cell_type = db$cellName[i],
                gene = c(p1, p2),
                direction = c(
                    rep("positive", length(p1)),
                    rep("negative", length(p2))
                ),
                level = db$level[i],
                stringsAsFactors = FALSE
            )
        })
        return(do.call(rbind, parts))
    }
    db
}

# ============================================================================
# Marker-finding methods
# ============================================================================

#' Find markers by fold change
#' @keywords internal
find_markers_fc <- function(
    exprs,
    clusters,
    top,
    min_log2fc,
    min_pct,
    pos_only,
    refs,
    max_pct_out
) {
    cts <- as.character(unique(clusters))
    cells <- names(clusters)
    markers1 <- setNames(rep("", length(cts)), cts)
    markers2 <- setNames(rep("", length(cts)), cts)
    n_dropped <- setNames(integer(length(cts)), cts)
    top_pos <- top[1]
    top_neg <- top[2]
    for (ct in cts) {
        ref_types <- if (is.null(refs)) NULL else refs[[ct]]
        if (!is.null(refs) && length(ref_types) == 0) {
            next  # `against` type left without other listed types
        }
        cells_in <- cells[clusters == ct]
        cells_out <- setdiff(cells, cells_in)
        m_in <- exprs[, cells_in, drop = FALSE]
        mean_in <- Matrix::rowMeans(m_in)
        pct_in <- Matrix::rowMeans(m_in > 0)
        rm(m_in)
        # one-vs-rest stats: the default reference group, or the pct_out
        # backing the max_pct_out guard / the negative candidates when
        # `against` is set
        m_out <- exprs[, cells_out, drop = FALSE]
        mean_out <- Matrix::rowMeans(m_out)
        pct_out <- Matrix::rowMeans(m_out > 0)
        rm(m_out)
        if (is.null(ref_types)) {
            log2fc <- log2((mean_in + 1e-6) / (mean_out + 1e-6))
            pct_cmp <- pct_out
        } else {
            # `against`: compute the stats against the reference group
            m_ref <- exprs[, cells[clusters %in% ref_types], drop = FALSE]
            mean_ref <- Matrix::rowMeans(m_ref)
            pct_ref <- Matrix::rowMeans(m_ref > 0)
            rm(m_ref)
            log2fc <- log2((mean_in + 1e-6) / (mean_ref + 1e-6))
            pct_cmp <- pct_ref
        }
        score <- log2fc * (pct_in - pct_cmp)

        keep <- pct_in >= min_pct
        if (pos_only) {
            keep <- keep & log2fc >= min_log2fc
        }
        # max_pct_out guard (positive candidates only): drop genes
        # expressed in most of ALL other cells, even when `against` is
        # set (the one-vs-rest pct_out is always computed)
        guard <- keep & pct_out > max_pct_out
        if (any(guard)) {
            keep[guard] <- FALSE
            n_dropped[[ct]] <- sum(guard)
        }
        idx <- which(keep)
        if (length(idx) > 0) {
            ord <- order(score[idx], decreasing = TRUE)
            take <- idx[ord[seq_len(min(top_pos, length(idx)))]]
            markers1[[ct]] <- paste(rownames(exprs)[take], collapse = ",")
        }

        if (!pos_only) {
            # negative markers: low in the type vs its comparison group
            idx2 <- which(log2fc <= -min_log2fc & pct_cmp >= min_pct)
            if (!is.null(ref_types)) {
                # with `against`: candidates must also be low in the type
                # vs the rest of the cells as in the default mode, so
                # that only sibling-specific negatives are kept
                log2fc_rest <- log2((mean_in + 1e-6) / (mean_out + 1e-6))
                idx2 <- intersect(
                    idx2,
                    which(log2fc_rest <= -min_log2fc & pct_out >= min_pct)
                )
            }
            if (length(idx2) > 0) {
                ord2 <- order(score[idx2], decreasing = FALSE)
                take2 <- idx2[ord2[seq_len(min(top_neg, length(idx2)))]]
                markers2[[ct]] <- paste(
                    rownames(exprs)[take2], collapse = ","
                )
            }
        }
    }
    if (sum(n_dropped) > 0) {
        dropped <- n_dropped[n_dropped > 0]
        warning(
            "Dropped ", sum(n_dropped), " positive-marker candidate(s) ",
            "expressed in more than ", max_pct_out * 100, "% of other ",
            "cells (max_pct_out): ",
            paste0(names(dropped), ": ", dropped, collapse = ", ")
        )
    }
    list(markers1 = markers1, markers2 = markers2)
}

#' Find markers with Seurat FindAllMarkers / FindMarkers
#' @keywords internal
find_markers_seurat <- function(
    exprs,
    clusters,
    top,
    min_log2fc,
    min_pct,
    pos_only,
    refs,
    max_pct_out
) {
    if (!inherits(exprs, "Seurat")) {
        stop(
            "method = \"seurat\" requires a Seurat object. ",
            "Please pass a Seurat object, or use method = ",
            "\"fc\" or \"presto\"."
        )
    }
    Seurat::Idents(exprs) <- unname(clusters[colnames(exprs)])
    cts <- as.character(unique(clusters))
    markers1 <- setNames(rep("", length(cts)), cts)
    markers2 <- setNames(rep("", length(cts)), cts)
    n_dropped <- setNames(integer(length(cts)), cts)
    top_pos <- top[1]
    top_neg <- top[2]

    # Guarded positive markers + negative markers for one type, given its
    # per-gene log2fc and its one-vs-rest pct_out
    extract_type_markers <- function(ct, genes, logfc, pct_out) {
        pos <- logfc >= 0
        if (any(pos)) {
            keep <- pct_out[pos] <= max_pct_out
            n_dropped[ct] <<- sum(!keep)
            cand <- which(pos)[keep]
            if (length(cand) > 0) {
                ord <- order(logfc[cand], decreasing = TRUE)
                take <- genes[
                    cand[ord[seq_len(min(top_pos, length(cand)))]]
                ]
                markers1[ct] <<- paste(take, collapse = ",")
            }
        }
        if (!pos_only) {
            neg <- which(logfc < 0)
            if (length(neg) > 0) {
                ord2 <- order(logfc[neg], decreasing = FALSE)
                take2 <- genes[neg[ord2[seq_len(min(top_neg, length(neg)))]]]
                markers2[ct] <<- paste(take2, collapse = ",")
            }
        }
    }

    if (is.null(refs)) {
        fam <- Seurat::FindAllMarkers(
            exprs,
            only.pos = pos_only,
            logfc.threshold = min_log2fc,
            min.pct = min_pct
        )
        fam <- fam[fam$gene %in% rownames(exprs), , drop = FALSE]
        if (nrow(fam) == 0) {
            stop("No marker genes found by Seurat::FindAllMarkers.")
        }
        logfc_col <- if ("avg_log2FC" %in% colnames(fam)) {
            "avg_log2FC"
        } else {
            "avg_logFC"
        }
        fam$cluster <- as.character(fam$cluster)
        for (ct in cts) {
            sub <- fam[fam$cluster == ct, , drop = FALSE]
            if (nrow(sub) > 0) {
                # pct.2 is the one-vs-rest pct_out (fraction of all the
                # other cells expressing the gene)
                extract_type_markers(
                    ct, sub$gene, sub[[logfc_col]], sub$pct.2
                )
            }
        }
    } else {
        # `against`: per-type two-group tests vs the reference group
        data_mat <- as.matrix_or_sparse(Seurat::GetAssayData(
            exprs, layer = "data"
        ))
        for (ct in cts) {
            ref_types <- refs[[ct]]
            if (length(ref_types) == 0) {
                next
            }
            fm <- Seurat::FindMarkers(
                exprs,
                ident.1 = ct,
                ident.2 = ref_types,
                only.pos = pos_only,
                logfc.threshold = min_log2fc,
                min.pct = min_pct
            )
            fm$gene <- rownames(fm)
            fm <- fm[fm$gene %in% rownames(exprs), , drop = FALSE]
            if (nrow(fm) == 0) {
                next
            }
            logfc_col <- if ("avg_log2FC" %in% colnames(fm)) {
                "avg_log2FC"
            } else {
                "avg_logFC"
            }
            # max_pct_out guard: the one-vs-rest pct_out over all cells
            # not of the type (computed even though the comparison is
            # against the reference group)
            cells_out <- colnames(data_mat)[clusters != ct]
            pct_out <- Matrix::rowMeans(
                data_mat[fm$gene, cells_out, drop = FALSE] > 0
            )
            extract_type_markers(ct, fm$gene, fm[[logfc_col]], pct_out)
        }
    }
    if (sum(n_dropped) > 0) {
        dropped <- n_dropped[n_dropped > 0]
        warning(
            "Dropped ", sum(n_dropped), " positive-marker candidate(s) ",
            "expressed in more than ", max_pct_out * 100, "% of other ",
            "cells (max_pct_out): ",
            paste0(names(dropped), ": ", dropped, collapse = ", ")
        )
    }
    list(markers1 = markers1, markers2 = markers2)
}

#' Find markers with presto wilcoxauc
#' @keywords internal
find_markers_presto <- function(
    exprs,
    clusters,
    top,
    min_log2fc,
    min_pct,
    pos_only,
    max_pct_out
) {
    if (!requireNamespace("presto", quietly = TRUE)) {
        stop(
            "Package 'presto' is required for method = 'presto'. ",
            "Install with: install.packages('presto')"
        )
    }
    # wilcoxauc's second positional argument is the group label vector
    # (`groups` in presto < 1.0, `y` in presto >= 1.0)
    res <- presto::wilcoxauc(as.matrix_or_sparse(exprs), clusters)
    res <- res[res$feature %in% rownames(exprs), , drop = FALSE]
    if (nrow(res) == 0) {
        stop("No marker genes found in the expression matrix rownames.")
    }
    cts <- as.character(unique(clusters))
    markers1 <- setNames(rep("", length(cts)), cts)
    markers2 <- setNames(rep("", length(cts)), cts)
    n_dropped <- setNames(integer(length(cts)), cts)
    top_pos <- top[1]
    top_neg <- top[2]
    res$group <- as.character(res$group)
    for (ct in cts) {
        sub <- res[res$group == ct, , drop = FALSE]
        if (nrow(sub) == 0) {
            next
        }
        # presto (>= 0.7) provides pct_in/pct_out; compute them from
        # exprs otherwise, as find_markers_fc does
        if (all(c("pct_in", "pct_out") %in% colnames(sub))) {
            # presto reports pct_in/pct_out as percentages (0-100);
            # rescale to fractions so the min_pct filter and the
            # logFC * (pct_in - pct_out) score match the fc backend
            pct_in <- sub$pct_in / 100
            pct_out <- sub$pct_out / 100
        } else {
            cells_in <- names(clusters)[clusters == ct]
            cells_out <- names(clusters)[clusters != ct]
            m_in <- exprs[sub$feature, cells_in, drop = FALSE]
            pct_in <- Matrix::rowMeans(m_in > 0)
            m_out <- exprs[sub$feature, cells_out, drop = FALSE]
            pct_out <- Matrix::rowMeans(m_out > 0)
        }
        # Same ranking as find_markers_fc: logFC * (pct_in - pct_out)
        score <- sub$logFC * (pct_in - pct_out)
        keep <- pct_in >= min_pct
        if (pos_only) {
            keep <- keep & sub$logFC >= min_log2fc
        }
        # max_pct_out guard on the positive candidates (one-vs-rest
        # pct_out)
        guard <- keep & pct_out > max_pct_out
        if (any(guard)) {
            keep[guard] <- FALSE
            n_dropped[[ct]] <- sum(guard)
        }
        idx <- which(keep)
        if (length(idx) > 0) {
            ord <- order(score[idx], decreasing = TRUE)
            take <- sub$feature[
                idx[ord[seq_len(min(top_pos, length(idx)))]]
            ]
            markers1[[ct]] <- paste(take, collapse = ",")
        }
        if (!pos_only) {
            idx2 <- which(sub$logFC <= -min_log2fc)
            if ("pct_out" %in% colnames(sub)) {
                idx2 <- idx2[pct_out[idx2] >= min_pct]
            }
            if (length(idx2) > 0) {
                ord2 <- order(sub$logFC[idx2], decreasing = FALSE)
                take2 <- sub$feature[
                    idx2[ord2[seq_len(min(top_neg, length(idx2)))]]
                ]
                markers2[[ct]] <- paste(take2, collapse = ",")
            }
        }
    }
    if (sum(n_dropped) > 0) {
        dropped <- n_dropped[n_dropped > 0]
        warning(
            "Dropped ", sum(n_dropped), " positive-marker candidate(s) ",
            "expressed in more than ", max_pct_out * 100, "% of other ",
            "cells (max_pct_out): ",
            paste0(names(dropped), ": ", dropped, collapse = ", ")
        )
    }
    list(markers1 = markers1, markers2 = markers2)
}

#' Coerce expression data to a matrix or a sparse matrix
#' @keywords internal
as.matrix_or_sparse <- function(x) {
    if (inherits(x, "dgTMatrix")) {
        return(methods::as(x, "dgCMatrix"))
    }
    x
}
