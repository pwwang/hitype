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
#' @param top Number of markers to return per cell type.
#' @param min_log2fc Minimum log2 fold change for a gene to be kept as a
#'  marker (when `pos_only = TRUE`).
#' @param min_pct Minimum fraction of cells in the cell type expressing
#'  the gene.
#' @param pos_only Only keep genes that are higher in the cell type than
#'  in the rest of the cells (positive markers). If `FALSE`, all genes
#'  passing `min_pct` are considered, ranked by score.
#' @param include_negative If `TRUE`, also fill the `geneSymbolmore2`
#'  column with the top down-regulated markers per cell type.
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
    top = 20,
    min_log2fc = 0.25,
    min_pct = 0.1,
    pos_only = TRUE,
    include_negative = FALSE,
    level = 1,
    format = c("universal", "db")
) {
    method <- match.arg(method)
    format <- match.arg(format)
    if (!is.numeric(top) || length(top) != 1 || is.na(top) || top < 1) {
        stop("`top` must be a single number >= 1")
    }
    if (!is.numeric(min_log2fc) || length(min_log2fc) != 1 ||
        is.na(min_log2fc) || min_log2fc < 0) {
        stop("`min_log2fc` must be a single number >= 0")
    }
    if (!is.numeric(min_pct) || length(min_pct) != 1 ||
        is.na(min_pct) || min_pct < 0 || min_pct > 1) {
        stop("`min_pct` must be a single number in [0, 1]")
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

    res <- switch(
        method,
        fc = find_markers_fc(
            exprs, clusters, top, min_log2fc, min_pct,
            pos_only, include_negative
        ),
        seurat = find_markers_seurat(
            exprs, clusters, top, min_log2fc, min_pct,
            pos_only, include_negative
        ),
        presto = find_markers_presto(
            exprs, clusters, top, min_log2fc, min_pct,
            pos_only, include_negative
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
    include_negative
) {
    cts <- as.character(unique(clusters))
    cells <- names(clusters)
    markers1 <- setNames(rep("", length(cts)), cts)
    markers2 <- setNames(rep("", length(cts)), cts)
    for (ct in cts) {
        cells_in <- cells[clusters == ct]
        cells_out <- setdiff(cells, cells_in)
        m_in <- exprs[, cells_in, drop = FALSE]
        mean_in <- Matrix::rowMeans(m_in)
        pct_in <- Matrix::rowMeans(m_in > 0)
        rm(m_in)
        m_out <- exprs[, cells_out, drop = FALSE]
        mean_out <- Matrix::rowMeans(m_out)
        pct_out <- Matrix::rowMeans(m_out > 0)
        rm(m_out)
        log2fc <- log2((mean_in + 1e-6) / (mean_out + 1e-6))
        score <- log2fc * (pct_in - pct_out)

        keep <- pct_in >= min_pct
        if (pos_only) {
            keep <- keep & log2fc >= min_log2fc
        }
        idx <- which(keep)
        if (length(idx) > 0) {
            ord <- order(score[idx], decreasing = TRUE)
            take <- idx[ord[seq_len(min(as.integer(top), length(idx)))]]
            markers1[[ct]] <- paste(rownames(exprs)[take], collapse = ",")
        }

        if (include_negative) {
            idx2 <- which(log2fc <= -min_log2fc & pct_out >= min_pct)
            if (length(idx2) > 0) {
                ord2 <- order(score[idx2], decreasing = FALSE)
                take2 <- idx2[ord2[seq_len(min(as.integer(top), length(idx2)))]]
                markers2[[ct]] <- paste(
                    rownames(exprs)[take2], collapse = ","
                )
            }
        }
    }
    list(markers1 = markers1, markers2 = markers2)
}

#' Find markers with Seurat FindAllMarkers
#' @keywords internal
find_markers_seurat <- function(
    exprs,
    clusters,
    top,
    min_log2fc,
    min_pct,
    pos_only,
    include_negative
) {
    if (!inherits(exprs, "Seurat")) {
        stop(
            "method = \"seurat\" requires a Seurat object. ",
            "Please pass a Seurat object, or use method = ",
            "\"fc\" or \"presto\"."
        )
    }
    Seurat::Idents(exprs) <- unname(clusters[colnames(exprs)])
    fam <- Seurat::FindAllMarkers(
        exprs,
        only.pos = pos_only || include_negative,
        logfc.threshold = min_log2fc,
        min.pct = min_pct
    )
    fam <- fam[fam$gene %in% rownames(exprs), , drop = FALSE]
    if (nrow(fam) == 0) {
        stop("No marker genes found by Seurat::FindAllMarkers.")
    }
    cts <- as.character(unique(clusters))
    markers1 <- setNames(rep("", length(cts)), cts)
    markers2 <- setNames(rep("", length(cts)), cts)
    logfc_col <- if ("avg_log2FC" %in% colnames(fam)) {
        "avg_log2FC"
    } else {
        "avg_logFC"
    }
    fam$cluster <- as.character(fam$cluster)
    for (ct in cts) {
        sub <- fam[fam$cluster == ct, , drop = FALSE]
        if (nrow(sub) > 0) {
            ord <- order(sub[[logfc_col]], decreasing = TRUE)
            take <- sub$gene[ord[seq_len(min(as.integer(top), nrow(sub)))]]
            markers1[[ct]] <- paste(take, collapse = ",")
        }
        if (include_negative && nrow(sub) > 0) {
            neg <- sub[sub[[logfc_col]] < 0, , drop = FALSE]
            if (nrow(neg) > 0) {
                ord2 <- order(neg[[logfc_col]], decreasing = FALSE)
                take2 <- neg$gene[
                    ord2[seq_len(min(as.integer(top), nrow(neg)))]
                ]
                markers2[[ct]] <- paste(take2, collapse = ",")
            }
        }
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
    include_negative
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
        idx <- which(keep)
        if (length(idx) > 0) {
            ord <- order(score[idx], decreasing = TRUE)
            take <- sub$feature[
                idx[ord[seq_len(min(as.integer(top), length(idx)))]]
            ]
            markers1[[ct]] <- paste(take, collapse = ",")
        }
        if (include_negative) {
            idx2 <- which(sub$logFC <= -min_log2fc)
            if ("pct_out" %in% colnames(sub)) {
                idx2 <- idx2[pct_out[idx2] >= min_pct]
            }
            if (length(idx2) > 0) {
                ord2 <- order(sub$logFC[idx2], decreasing = FALSE)
                take2 <- sub$feature[
                    idx2[ord2[seq_len(min(as.integer(top), length(idx2)))]]
                ]
                markers2[[ct]] <- paste(take2, collapse = ",")
            }
        }
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
