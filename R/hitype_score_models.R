#' Score cells with the models trained by train_weights()
#'
#' Scores each cell with the linear predictors of the per-cell-type models
#' persisted by [train_weights()] with `return_models = TRUE` (methods
#' `"glmnet"` and `"lr"` only): for a cell with expression `x` of the model
#' features, the score of cell type `t` is
#' `eta_t = (Intercept)_t + sum_f coef_tf * (x_f - center_f) / scale_f`,
#' where `center`/`scale` are the per-gene centering/scaling recorded when
#' the models were trained. Genes of the model features missing from
#' `exprs` contribute 0; extra genes of `exprs` are ignored.
#'
#' The centering/scaling is folded into the coefficients, so a sparse
#' `exprs` matrix is never densified by the `(x - center) / scale` shift:
#' `eta_t = (Intercept)_t - sum_f b_tf * center_f + sum_f b_tf * x_f` with
#' `b_tf = coef_tf / scale_f`.
#'
#' @param exprs Input scRNA-seq expression matrix (genes x cells, the same
#'  convention as [hitype_score()]).
#' @param models The model bundle returned by
#'  [train_weights()] with `return_models = TRUE` (the `models` element).
#' @param margin Cells whose top-minus-second score (`margins`) is below
#'  `margin` are assigned `"Unknown"` instead of their top cell type.
#'  `0` (default) assigns every cell.
#'
#' @return A list with:
#'  \describe{
#'    \item{`scores`}{A matrix (cells x cell types) of the linear
#'      predictors, with the cell types in the order of the model bundle.}
#'    \item{`assignments`}{A named vector with the top-scoring cell type
#'      of every cell (or `"Unknown"` for cells below the `margin`).}
#'    \item{`margins`}{A named numeric vector with the top-minus-second
#'      score of every cell.}
#'  }
#' @export
hitype_score_models <- function(exprs, models, margin = 0) {
    if (is.data.frame(exprs)) {
        exprs <- as.matrix(exprs)
    }
    if (
        !is.matrix(exprs) &&
        !any(class(exprs) %in% c("dgCMatrix", "dgTMatrix"))
    ) {
        stop(
            "Input scRNA-seq data must be a matrix or data.frame ",
            "(genes x cells)"
        )
    }
    if (sum(dim(exprs)) == 0) {
        stop("Input scRNA-seq data is empty")
    }
    if (!is.numeric(margin) || length(margin) != 1) {
        stop("`margin` must be a single number")
    }
    if (!is.list(models) || is.null(models$coefs)) {
        stop("model-based scoring is only available for methods glmnet and lr")
    }
    coefs <- models$coefs
    if (!is.list(coefs) || length(coefs) == 0) {
        stop("The model bundle contains no cell type models")
    }

    if (isTRUE(models$scaled_input)) {
        warning(
            "The models were trained on already-scaled input ",
            "(`scaled = TRUE`); `exprs` is expected to be scaled too, ",
            "so no centering/scaling is applied.",
            immediate. = TRUE
        )
    }

    features <- models$features
    present <- intersect(features, rownames(exprs))
    if (length(present) == 0) {
        stop("No model features are present in `exprs`")
    }
    # Genes missing from `exprs` contribute 0: keep only the present rows,
    # the per-type dot product below simply does not include the others
    z <- exprs[present, , drop = FALSE]
    cells <- colnames(z)
    ncell <- ncol(z)

    if (isTRUE(models$scaled_input)) {
        # Already-scaled input: b = coef and the center correction is 0
        center <- setNames(rep(0, length(present)), present)
        scale_v <- setNames(rep(1, length(present)), present)
    } else {
        center <- models$center[present]
        scale_v <- models$scale[present]
    }

    # Per-type coefficient row over the present features, centering and
    # scaling folded in (b_tf = coef_tf / scale_f) so the multiplication
    # below stays sparse-friendly:
    #   sum_f coef_tf * (x_f - center_f) / scale_f
    # = sum_f b_tf * x_f - sum_f b_tf * center_f
    types <- names(coefs)
    B <- do.call(rbind, lapply(coefs, function(cf) {
        b <- cf[present] / scale_v
        b[is.na(b)] <- 0  # a feature absent from the fit coefficients
        b
    }))
    intercepts <- vapply(coefs, function(cf) {
        b0 <- cf["(Intercept)"]
        if (is.na(b0)) 0 else unname(b0)
    }, numeric(1))

    scores <- as.matrix(t(z) %*% t(B))  # cells x types
    colnames(scores) <- types
    # per type: eta = intercept - center correction + data term
    scores <- sweep(scores, 2, intercepts - as.numeric(B %*% center), "+")
    if (!is.null(cells)) {
        rownames(scores) <- cells
    }

    ord <- as.matrix(apply(scores, 1, order, decreasing = TRUE))
    best1 <- scores[cbind(seq_len(ncell), ord[1, ])]
    if (ncol(scores) == 1) {
        # A single cell type: no runner-up, margins are infinite
        best2 <- rep(-Inf, ncell)
    } else {
        best2 <- scores[cbind(seq_len(ncell), ord[2, ])]
    }
    assignments <- colnames(scores)[ord[1, ]]
    margins <- best1 - best2
    if (!is.null(cells)) {
        names(assignments) <- cells
        names(margins) <- cells
    }
    assignments[margins < margin] <- "Unknown"

    list(scores = scores, assignments = assignments, margins = margins)
}
