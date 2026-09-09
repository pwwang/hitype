# Model-based scoring: train_weights(return_models = TRUE) persists the
# per-type linear models (glmnet/lr) and hitype_score_models() scores
# cells with the model equations (linear predictors), with a margin-based
# Unknown threshold. The default return_models = FALSE return value is
# unchanged.

skip_if_not_installed("glmnet")

# Synthetic fixture: 3 well-separated cell types (T1/T2/T3) of 3 marker
# genes each over a noisy background. The training and held-out cells come
# from the same generation process but the held-out cells never enter
# train_weights().
synth_fixture <- function(seed = 123, train_per = 30, holdout_per = 10) {
    types <- paste0("T", 1:3)
    marker_genes <- list(
        T1 = paste0("G", 1:3), T2 = paste0("G", 4:6), T3 = paste0("G", 7:9)
    )
    per <- train_per + holdout_per
    set.seed(seed)
    exprs_all <- matrix(
        rnorm(9 * 3 * per, 0, 0.3), nrow = 9,
        dimnames = list(paste0("G", 1:9), paste0("c", 1:(3 * per)))
    )
    for (i in seq_along(types)) {
        cells <- ((i - 1) * per + 1):(i * per)
        exprs_all[marker_genes[[i]], cells] <-
            exprs_all[marker_genes[[i]], cells] + 4
    }
    train_cells <- unlist(lapply(seq_along(types), function(i) {
        ((i - 1) * per + 1):((i - 1) * per + train_per)
    }))
    hold_cells <- setdiff(seq_len(ncol(exprs_all)), train_cells)
    clusters <- rep(types, each = train_per)
    names(clusters) <- colnames(exprs_all)[train_cells]
    list(
        file_markers = data.frame(
            cell_type = rep(types, each = 3),
            gene = paste0("G", 1:9),
            stringsAsFactors = FALSE
        ),
        exprs = exprs_all[, train_cells, drop = FALSE],
        clusters = clusters,
        holdout = exprs_all[, hold_cells, drop = FALSE],
        truth = rep(types, each = holdout_per)
    )
}

test_that("glmnet model bundle scores held-out cells correctly", {
    fx <- synth_fixture()
    res <- train_weights(
        path_to_gs = fx$file_markers,
        exprs = fx$exprs,
        clusters = fx$clusters,
        method = "glmnet",
        data_split = c(1),  # the fixture holds the test cells out itself
        return_models = TRUE
    )
    expect_equal(names(res), c("weights", "models", "labels"))
    expect_equal(res$labels, c("T1", "T2", "T3"))
    m <- res$models
    expect_equal(m$method, "glmnet")
    expect_equal(m$level, 1L)
    expect_equal(m$features, paste0("G", 1:9))
    expect_equal(names(m$center), m$features)
    expect_equal(names(m$scale), m$features)
    expect_null(m$scaled_input)  # scaled = FALSE: no scaled-input flag
    expect_equal(names(m$coefs), res$labels)
    for (cf in m$coefs) {
        expect_equal(names(cf), c("(Intercept)", m$features))
    }

    sc <- hitype_score_models(fx$holdout, m)
    expect_equal(dim(sc$scores), c(30, 3))
    expect_equal(colnames(sc$scores), res$labels)
    expect_equal(names(sc$assignments), colnames(fx$holdout))
    # Well-separated types: allow a couple of errors at most, and every
    # confidently assigned cell has a positive top-minus-second margin
    expect_true(sum(sc$assignments == fx$truth) >= 28)
    expect_true(all(sc$margins > 0))
})

test_that("margin threshold sends uncertain cells to Unknown", {
    fx <- synth_fixture()
    res <- train_weights(
        path_to_gs = fx$file_markers,
        exprs = fx$exprs,
        clusters = fx$clusters,
        method = "glmnet",
        data_split = c(1),
        return_models = TRUE
    )
    sc0 <- hitype_score_models(fx$holdout, res$models, margin = 0)
    sc_all <- hitype_score_models(fx$holdout, res$models, margin = 1e6)
    expect_equal(sum(sc0$assignments == "Unknown"), 0)
    expect_equal(sum(sc_all$assignments == "Unknown"), length(fx$truth))
})

test_that("a feature missing from exprs is tolerated", {
    fx <- synth_fixture()
    res <- train_weights(
        path_to_gs = fx$file_markers,
        exprs = fx$exprs,
        clusters = fx$clusters,
        method = "glmnet",
        data_split = c(1),
        return_models = TRUE
    )
    sc_full <- hitype_score_models(fx$holdout, res$models)
    # Drop G1 (a T1 marker): the function runs and the assignments stay
    # nearly unchanged
    sc_missing <- hitype_score_models(
        fx$holdout[setdiff(rownames(fx$holdout), "G1"), , drop = FALSE],
        res$models
    )
    expect_equal(length(sc_missing$assignments), length(fx$truth))
    expect_true(sum(sc_missing$assignments == sc_full$assignments) >= 27)
})

test_that("lr bundle has the same shape and scores", {
    fx <- synth_fixture()
    res <- train_weights(
        path_to_gs = fx$file_markers,
        exprs = fx$exprs,
        clusters = fx$clusters,
        method = "lr",
        data_split = c(1),
        return_models = TRUE
    )
    expect_equal(names(res), c("weights", "models", "labels"))
    expect_equal(res$labels, c("T1", "T2", "T3"))
    expect_true(all(vapply(res$models$coefs, function(cf) {
        "(Intercept)" %in% names(cf)
    }, logical(1))))
    sc <- hitype_score_models(fx$holdout, res$models)
    expect_equal(dim(sc$scores), c(30, 3))
    expect_equal(colnames(sc$scores), res$labels)
    expect_equal(names(sc$assignments), colnames(fx$holdout))
    expect_equal(length(sc$margins), 30)
    expect_true(sum(sc$assignments == fx$truth) >= 28)
})

test_that("return_models = FALSE keeps the historical return value", {
    fx <- synth_fixture()
    args <- list(
        path_to_gs = fx$file_markers,
        exprs = fx$exprs,
        clusters = fx$clusters,
        method = "glmnet",
        data_split = c(1)
    )
    w <- do.call(train_weights, args)
    w2 <- do.call(train_weights, c(args, list(return_models = FALSE)))
    expect_identical(w, w2)
    # Historical structure: universal-format weights data.frame
    expect_s3_class(w, "data.frame")
    expect_equal(
        colnames(w),
        c("cell_type", "gene", "direction", "weight", "level")
    )
    expect_equal(w$level, rep(1L, nrow(w)))
    expect_true(all(w$direction %in% c("positive", "negative")))
    expect_setequal(unique(w$cell_type), c("T1", "T2", "T3"))
})

test_that("unsupported methods error at scoring time", {
    fx <- synth_fixture()
    res <- train_weights(
        path_to_gs = fx$file_markers,
        exprs = fx$exprs,
        clusters = fx$clusters,
        method = "uniform",
        data_split = c(1),
        return_models = TRUE
    )
    expect_null(res$models$coefs)
    expect_null(res$labels)
    expect_equal(
        colnames(res$weights),
        c("cell_type", "gene", "direction", "weight", "level")
    )
    expect_error(
        hitype_score_models(fx$holdout, res$models),
        "only available for methods glmnet and lr"
    )
})
