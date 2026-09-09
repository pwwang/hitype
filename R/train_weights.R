#' Train weights for the markers
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom stats median
#' @importFrom stats na.omit
#'
#' @param path_to_gs Path to the gene set file without weights. The
#'  training cell types are the cell types of the cells (see `clusters`):
#'  a cell type in the marker file with cells of that type in the data is
#'  trained on its own markers in the file. A cell type in the marker file
#'  with no cells of that type in the data is ignored with a warning, and
#'  its markers are pooled for the cell types in the data that are not
#'  covered by the marker file (each of them is trained on the pooled
#'  markers). Cell types in the data that are neither covered by the marker
#'  file nor by pooled markers are not trained (with a warning).
#' @param exprs The expression matrix, or a seurat object
#'  (rows: genes, columns: samples/cells)
#' @param level The level of the gene sets to train weights for if
#'  you have multiple levels of gene sets.
#' @param scaled Whether the expression matrix is scaled
#' @param clusters The cell types of the cells. When `exprs` is a Seurat
#'  object, it can be `NULL` (default) to take the cell types from
#'  `Seurat::Idents()`, or a column name in the `meta.data` of the Seurat
#'  object that holds the cell type of each cell. When `exprs` is a matrix,
#'  it should be a named vector of cell types (names - cell names).
#' @param range The quantization range for the db output only (`format =
#'  "db"`): the 4-element form `c(-low, -high, low, high)` (default
#'  `c(-5, -1, 1, 5)`) scales the positive weights into `[low, high]` and
#'  the negative ones into `[-high, -low]` before rounding them to
#'  integers, so genuinely positive markers always stay positive. The
#'  universal output (default) keeps the raw trained weights as-is.
#' @param data_split A vector of fractions for training, validation and
#'  testing. If only two fractions are provided, no testing set will be used.
#' @param epochs The number of epochs to train (lrp method only)
#' @param batch_size The batch size (lrp method only)
#' @param run_weights_on_test Whether to run the weights on the test set.
#'  Requires that `data_split` has three elements.
#' @param cv_folds Number of cross-validation folds for weight estimation.
#'  When > 1, weights are averaged across folds for stability.
#'  Default is 1 (no cross-validation). Used by lr, glmnet, and lrp methods.
#' @param method The weight learning method. One of:
#'  \describe{
#'    \item{"uniform"}{All markers get equal weight (= 1). Fast baseline.}
#'    \item{"correlation"}{Pearson correlation between each marker and
#'      the binary cluster indicator. Simple, interpretable.}
#'    \item{"lr"}{Logistic regression coefficients (one-vs-rest).
#'      Classic ML approach. Use cv_folds for stability.}
#'    \item{"glmnet"}{Sparse logistic regression with elastic net penalty
#'      (alpha = 0.5). Automatically zeros out uninformative markers.
#'      Recommended for most users. Requires the glmnet package.}
#'    \item{"rf"}{Random forest permutation importance. Captures non-linear
#'      marker interactions. Requires the ranger package.}
#'    \item{"xgb"}{XGBoost gain-based feature importance.
#'      State-of-the-art tree method. Requires the xgboost package.}
#'    \item{"lrp"}{Neural network + Layer-wise Relevance Propagation.
#'      Deep learning approach. Requires the keras and innsight packages.}
#'  }
#' @param format The format of the output data frame. One of
#'   `"universal"` (default) or `"db"` (the hitype/ScType wide format).
#' @param drop_zero Whether to drop markers whose trained weight is exactly
#'   zero from the output (default `TRUE`). Such markers carry no learned
#'   direction — with the glmnet method (coefficients at `lambda.1se`) most
#'   candidate markers are zeroed. Dropping them makes zero mean "no
#'   direction": the marker is absent from the cell type's list. Pass
#'   `FALSE` to keep every candidate marker (with `format = "db"` they are
#'   then encoded as `*`).
#' @param pos_only Whether to keep only the markers with a positive
#'   trained weight in the output (default `FALSE`). A marker can be trained
#'   with a negative weight (anti-correlating with the cell type, e.g. by
#'   the `"correlation"` method); such rows are reported with
#'   `direction = "negative"` and are dropped when `pos_only = TRUE`,
#'   so the returned table contains markers overexpressed in each cell type
#'   only.
#' @param return_models Whether to also return the fitted per-cell-type
#'   prediction models (default `FALSE`). With `FALSE` (default) the return
#'   value is unchanged: a data frame with the weights. With `TRUE`, a list
#'   is returned instead:
#'   \describe{
#'     \item{`weights`}{The usual weights data frame (universal or db
#'       format).}
#'     \item{`models`}{A model bundle list that [hitype_score_models()]
#'       can score new cells with:
#'       \describe{
#'         \item{`method`}{The weight learning method.}
#'         \item{`level`}{The level of the gene sets.}
#'         \item{`features`}{The ordered gene list the models were fit on.}
#'         \item{`center`, `scale`}{Numeric vectors named by `features`
#'           with the per-gene centering/scaling that was applied to the
#'           training matrix, so scoring can reproduce it exactly.}
#'         \item{`coefs`}{A named list, one entry per trained cell type,
#'           each a named numeric coefficient vector of the linear
#'           predictor including the intercept as `"(Intercept)"`. `NULL`
#'           for methods without linear models.}
#'       }}
#'     \item{`labels`}{The trained cell types (names of `coefs`).}
#'   }
#'   Model bundles are only produced for the linear-model methods
#'   `"glmnet"` and `"lr"`; with any other method `coefs` is `NULL`,
#'   `labels` is `NULL`, and [hitype_score_models()] errors.
#' @param seed Random seed for reproducibility
#' @return By default, a data frame with the weights in the universal marker
#'  format (default) or the db format (`format = "db"`), that can be used
#'  directly by [gs_prepare()]. With `return_models = TRUE`, a list with the
#'  weights data frame, the model bundle and the trained cell type labels
#'  (see `return_models`).
#'
#' @export
train_weights <- function(
    path_to_gs,
    exprs,
    level = 1,
    scaled = FALSE,
    clusters = NULL,
    range = c(-5, -1, 1, 5),
    data_split = c(0.7, 0.2, 0.1),
    epochs = 20,
    batch_size = 32,
    run_weights_on_test = TRUE,
    cv_folds = 1,
    method = c("glmnet", "lr", "rf", "xgb", "lrp", "correlation", "uniform"),
    format = c("universal", "db"),
    drop_zero = TRUE,
    pos_only = FALSE,
    return_models = FALSE,
    seed = 8525
) {
    method <- match.arg(method)
    format <- match.arg(format)
    set.seed(seed)

    data <- prepare_data_for_training(
        path_to_gs, exprs, level, scaled, clusters
    )

    clusters <- as.character(data$clusters)
    # Keep the cluster types consistent: the weight-learning methods iterate
    # `unique(data$clusters)` and compare it against the `clusters` argument
    # (`clusters == ct`), and `compile_weights()` must match `output_node`
    # against the gene-set cell-type names. So both must carry the cell-type
    # labels as characters. When `exprs` is a Seurat object, `data$clusters`
    # is a factor: converting only the local to integer codes made
    # `integer == factor` all-FALSE — a silently all-zero response (every
    # `cv.glmnet` fit then failed and was swallowed by the tryCatch,
    # yielding the rescale midpoint for all weights).
    data$clusters <- clusters
    uclusters <- unique(clusters)

    # Hold-out test set (used only if data_split has 3 elements)
    test_data_x <- NULL
    test_data_y <- NULL
    if (length(data_split) == 3) {
        test_idx <- sample(
            seq_len(nrow(data$z)), floor(nrow(data$z) * data_split[3])
        )
        test_data_x <- data$z[test_idx, , drop = FALSE]
        test_data_y <- clusters[test_idx]
        rest_idx <- setdiff(seq_len(nrow(data$z)), test_idx)
        data$z <- data$z[rest_idx, , drop = FALSE]
        clusters <- clusters[rest_idx]
    }

    result <- switch(
        method,
        uniform    = train_uniform(data, clusters),
        correlation = train_correlation(data, clusters),
        lr         = train_lr(data, clusters, cv_folds, return_models),
        glmnet     = train_glmnet(data, clusters, cv_folds, return_models),
        rf         = train_rf(data, clusters),
        xgb        = train_xgb(data, clusters),
        lrp        = train_lrp(
            data, clusters, uclusters, cv_folds,
            epochs, batch_size, data_split
        )
    )

    # The linear-model methods (glmnet, lr) return the fitted per-type
    # coefficient vectors alongside the weights table when models are
    # requested, so the bundle and the table come from the same fits.
    model_coefs <- NULL
    if (return_models && method %in% c("glmnet", "lr")) {
        model_coefs <- result$coefs
        result <- result$weights
    } else if (return_models && method %in% c("rf", "xgb", "lrp")) {
        warning(
            paste0(
                "`return_models = TRUE` is not supported for method '",
                method, "'; model-based scoring is only available for ",
                "methods 'glmnet' and 'lr'"
            ),
            immediate. = TRUE
        )
    }

    weights <- compile_weights(
        result, data$gs, level, range, format, drop_zero, pos_only
    )
    if (!is.null(test_data_x) && run_weights_on_test) {
        run_weights_on_test_data(
            weights, exprs, clusters, scaled, rownames(test_data_x)
        )
    }
    if (!return_models) {
        return(weights)
    }

    # Report the cell types of the weights table (in its order): fits of
    # types the compiled weights drop (no markers left) are not included
    if (!is.null(model_coefs)) {
        model_coefs <- model_coefs[
            intersect(names(data$gs), names(model_coefs))
        ]
    }
    bundle <- list(
        method = method,
        level = as.integer(level),
        features = colnames(data$z),
        center = data$center,
        scale = data$scale,
        coefs = model_coefs
    )
    if (isTRUE(data$scaled_input)) bundle$scaled_input <- TRUE
    if (is.null(model_coefs)) bundle$models <- list()
    list(weights = weights, models = bundle, labels = names(model_coefs))
}

# ============================================================================
# Weight-learning methods
# ============================================================================

#' Uniform weights (baseline)
#' @keywords internal
train_uniform <- function(data, clusters) {
    uctypes <- unique(data$clusters)
    markers <- colnames(data$z)
    expand.grid(
        output_node = uctypes,
        feature = markers,
        value = 1,
        stringsAsFactors = FALSE
    )
}

#' Correlation-based weights
#' @keywords internal
#' @importFrom stats cor
train_correlation <- function(data, clusters) {
    uctypes <- unique(data$clusters)
    markers <- colnames(data$z)
    out <- lapply(uctypes, function(ct) {
        y <- as.numeric(clusters == ct)
        vals <- apply(data$z, 2, function(x) suppressWarnings(cor(x, y)))
        vals[is.na(vals)] <- 0
        data.frame(
            output_node = ct,
            feature = markers,
            value = vals,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

#' Logistic regression weights (one-vs-rest)
#' @keywords internal
#' @importFrom stats glm binomial coef
train_lr <- function(data, clusters, cv_folds = 1, return_models = FALSE) {
    uctypes <- unique(data$clusters)
    markers <- colnames(data$z)
    z <- data$z

    if (cv_folds > 1) {
        n <- nrow(z)
        fold <- integer(n)
        for (ct in clusters) {
            ct_idx <- which(clusters == ct)
            fold[ct_idx] <- sample(rep_len(seq_len(cv_folds), length(ct_idx)))
        }
        res <- lapply(uctypes, function(ct) {
            y <- as.numeric(clusters == ct)
            fold_fits <- lapply(seq_len(cv_folds), function(k) {
                train_idx <- which(fold != k)
                z_train <- z[train_idx, , drop = FALSE]
                y_train <- y[train_idx]
                fit <- suppressWarnings(
                    glm(y_train ~ ., data = as.data.frame(z_train),
                        family = binomial())
                )
                cf <- coef(fit)[-1]  # drop intercept
                cf[is.na(cf)] <- 0
                list(cf = cf, mcf = model_coefs_lr(fit, markers, return_models))
            })
            vals <- Reduce(`+`, lapply(fold_fits, `[[`, "cf")) / cv_folds
            list(
                df = data.frame(
                    output_node = ct,
                    feature = markers,
                    value = vals,
                    stringsAsFactors = FALSE
                ),
                mcoefs = fold_models_lr(fold_fits, return_models)
            )
        })
    } else {
        res <- lapply(uctypes, function(ct) {
            y <- as.numeric(clusters == ct)
            fit <- suppressWarnings(
                glm(y ~ ., data = as.data.frame(z), family = binomial())
            )
            vals <- coef(fit)[-1]
            vals[is.na(vals)] <- 0
            list(
                df = data.frame(
                    output_node = ct,
                    feature = markers,
                    value = vals,
                    stringsAsFactors = FALSE
                ),
                mcoefs = model_coefs_lr(fit, markers, return_models)
            )
        })
    }
    collect_models_result(res, uctypes, return_models)
}

#' Extract the per-type linear predictor coefficients of a glm fit
#'
#' Aligned by name to c("(Intercept)", markers): an aliased term can be
#' absent from `coef()`, it then stays 0.
#' @keywords internal
model_coefs_lr <- function(fit, markers, return_models) {
    if (!return_models) {
        return(NULL)
    }
    mcf <- coef(fit)
    mcf[is.na(mcf)] <- 0
    out <- setNames(
        c(0, rep(0, length(markers))), c("(Intercept)", markers)
    )
    out[names(mcf)] <- mcf
    out
}

#' Average per-fold linear predictor coefficients (cv_folds > 1)
#' @keywords internal
fold_models_lr <- function(fold_fits, return_models) {
    if (!return_models) {
        return(NULL)
    }
    Reduce(`+`, lapply(fold_fits, `[[`, "mcf")) / length(fold_fits)
}

#' Turn per-type fit results into the weights table and (optionally) the
#' coefficient list of the model bundle
#' @keywords internal
collect_models_result <- function(res, uctypes, return_models) {
    df <- do.call(rbind, lapply(res, `[[`, "df"))
    if (!return_models) {
        return(df)
    }
    coefs <- setNames(lapply(res, `[[`, "mcoefs"), uctypes)
    # A type whose fit failed (mcoefs NULL) is not in the bundle
    coefs <- coefs[!vapply(coefs, is.null, logical(1))]
    list(weights = df, coefs = coefs)
}

#' Sparse logistic regression weights (glmnet)
#' @keywords internal
train_glmnet <- function(data, clusters, cv_folds = 5, return_models = FALSE) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' is required for method='glmnet'. ",
             "Install with: install.packages('glmnet')")
    }
    uctypes <- unique(data$clusters)
    markers <- colnames(data$z)
    z <- data$z
    # The weights table takes the coefficients at lambda.1se; the model
    # bundle the coefficients (with the intercept) at lambda.min — both
    # from the same fit so they never diverge.
    glmnet_model_coefs <- function(fit) {
        if (!return_models) {
            return(NULL)
        }
        setNames(
            as.numeric(coef(fit, s = "lambda.min")),
            c("(Intercept)", markers)
        )
    }

    if (cv_folds > 1) {
        n <- nrow(z)
        fold <- integer(n)
        for (ct_val in clusters) {
            ct_idx <- which(clusters == ct_val)
            fold[ct_idx] <- sample(
                rep_len(seq_len(cv_folds), length(ct_idx))
            )
        }
        res <- lapply(uctypes, function(ct) {
            y <- as.numeric(clusters == ct)
            fold_fits <- lapply(seq_len(cv_folds), function(k) {
                train_idx <- which(fold != k)
                z_train <- as.matrix(z[train_idx, , drop = FALSE])
                y_train <- y[train_idx]
                fit <- tryCatch(
                    glmnet::cv.glmnet(
                        z_train, y_train,
                        family = "binomial", alpha = 0.5
                    ),
                    error = function(e) NULL
                )
                if (is.null(fit)) {
                    return(list(
                        cf = setNames(rep(0, length(markers)), markers),
                        mcf = if (return_models) {
                            setNames(
                                rep(0, length(markers) + 1),
                                c("(Intercept)", markers)
                            )
                        } else {
                            NULL
                        },
                        ok = FALSE
                    ))
                }
                lam <- if (inherits(fit, "cv.glmnet")) fit$lambda.1se
                       else median(fit$lambda)
                cf <- as.numeric(coef(fit, s = lam))[-1]
                cf[is.na(cf)] <- 0
                list(cf = cf, mcf = glmnet_model_coefs(fit), ok = TRUE)
            })
            vals <- Reduce(`+`, lapply(fold_fits, `[[`, "cf")) / cv_folds
            # A type whose fits all failed is left out of the bundle (the
            # weights table path handles it with all-zero coefficients)
            mcoefs <- NULL
            if (return_models) {
                ok <- vapply(fold_fits, function(ff) ff$ok, logical(1))
                if (any(ok)) {
                    mcoefs <- Reduce(
                        `+`, lapply(fold_fits, `[[`, "mcf")
                    ) / cv_folds
                }
            }
            list(
                df = data.frame(
                    output_node = ct,
                    feature = markers,
                    value = vals,
                    stringsAsFactors = FALSE
                ),
                mcoefs = mcoefs
            )
        })
        return(collect_models_result(res, uctypes, return_models))
    }

    res <- lapply(uctypes, function(ct) {
        y <- as.numeric(clusters == ct)
        zm <- as.matrix(z)
        fit <- tryCatch(
            glmnet::cv.glmnet(
                zm, y, family = "binomial", alpha = 0.5
            ),
            error = function(e) NULL
        )
        mcoefs <- NULL
        if (is.null(fit)) {
            vals <- rep(0, length(markers))
        } else {
            lam <- if (inherits(fit, "cv.glmnet")) fit$lambda.1se
                   else median(fit$lambda)
            vals <- as.numeric(coef(fit, s = lam))[-1]
            mcoefs <- glmnet_model_coefs(fit)
        }
        vals[is.na(vals)] <- 0
        list(
            df = data.frame(
                output_node = ct,
                feature = markers,
                value = vals,
                stringsAsFactors = FALSE
            ),
            mcoefs = mcoefs
        )
    })
    collect_models_result(res, uctypes, return_models)
}

#' Random forest permutation importance weights
#' @keywords internal
train_rf <- function(data, clusters) {
    if (!requireNamespace("ranger", quietly = TRUE)) {
        stop("Package 'ranger' is required for method='rf'. ",
             "Install with: install.packages('ranger')")
    }
    uctypes <- unique(data$clusters)
    markers <- colnames(data$z)
    zdf <- as.data.frame(data$z)
    zdf$cluster <- factor(clusters)

    fit <- ranger::ranger(
        cluster ~ ., data = zdf,
        importance = "permutation",
        num.trees = 500
    )
    imp <- ranger::importance(fit)
    imp[is.na(imp)] <- 0

    out <- lapply(uctypes, function(ct) {
        # Overall importance — same values for all cell types
        data.frame(
            output_node = ct,
            feature = markers,
            value = imp,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

#' XGBoost gain-based importance weights
#' @keywords internal
train_xgb <- function(data, clusters) {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop("Package 'xgboost' is required for method='xgb'. ",
             "Install with: install.packages('xgboost')")
    }
    uctypes <- unique(data$clusters)
    markers <- colnames(data$z)
    uclusters_int <- as.integer(factor(clusters)) - 1  # 0-based
    n_classes <- length(uctypes)

    dtrain <- xgboost::xgb.DMatrix(
        data = as.matrix(data$z),
        label = uclusters_int
    )
    params <- list(
        objective = "multi:softprob",
        num_class = n_classes,
        max_depth = 6,
        eta = 0.3,
        nthread = 1,
        verbosity = 0
    )
    fit <- xgboost::xgb.train(
        params = params, data = dtrain, nrounds = 50
    )
    imp <- xgboost::xgb.importance(
        feature_names = markers, model = fit
    )
    # Build complete matrix: all features × all cell types
    imp_vec <- setNames(rep(0, length(markers)), markers)
    imp_vec[imp$Feature] <- imp$Gain
    imp_vec[is.na(imp_vec)] <- 0

    out <- lapply(uctypes, function(ct) {
        data.frame(
            output_node = ct,
            feature = markers,
            value = imp_vec,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, out)
}

#' LRP-based weights (neural network)
#' @keywords internal
train_lrp <- function(
    data, clusters, uclusters, cv_folds,
    epochs, batch_size, data_split
) {
    if (!requireNamespace("keras", quietly = TRUE)) {
        stop("Package 'keras' is required for method='lrp'. ",
             "Install with: install.packages('keras'); ",
             "keras::install_keras()")
    }
    if (!requireNamespace("innsight", quietly = TRUE)) {
        stop("Package 'innsight' is required for method='lrp'. ",
             "Install with: install.packages('innsight')")
    }

    if (cv_folds > 1) {
        cv_clusters <- clusters
        n_cv <- nrow(data$z)
        fold <- integer(n_cv)
        for (ct in uclusters) {
            ct_idx <- which(cv_clusters == ct)
            fold[ct_idx] <- sample(rep_len(seq_len(cv_folds), length(ct_idx)))
        }
        all_results <- list()
        for (k in seq_len(cv_folds)) {
            cat(sprintf("CV fold %d/%d\n", k, cv_folds))
            train_idx <- which(fold != k)
            test_fold_idx <- which(fold == k)
            all_results[[k]] <- train_lrp_single(
                data$z, cv_clusters, uclusters,
                train_idx, test_fold_idx,
                epochs, batch_size, data_split
            )
        }
        result <- do.call(rbind, all_results)
        result <- result %>%
            group_by(output_node, feature) %>%
            summarise(value = mean(value), .groups = "drop")
    } else {
        result <- train_lrp_single(
            data$z, clusters, uclusters,
            seq_len(nrow(data$z)), seq_len(nrow(data$z)),
            epochs, batch_size, data_split
        )
    }
    result
}

#' Train a single LRP model (helper)
#' @keywords internal
train_lrp_single <- function(
    z, clusters, uclusters,
    train_idx, test_idx,
    epochs, batch_size, data_split
) {
    model <- keras::keras_model_sequential() %>%
        keras::layer_dense(
            units = 64, activation = "relu",
            input_shape = ncol(z)
        ) %>%
        keras::layer_dropout(rate = 0.2) %>%
        keras::layer_dense(units = 64, activation = "relu") %>%
        keras::layer_dropout(rate = 0.2) %>%
        keras::layer_dense(
            units = length(uclusters), activation = "softmax"
        )

    model %>% keras::compile(
        loss = "categorical_crossentropy",
        optimizer = "adam",
        metrics = c("accuracy")
    )
    model %>% keras::fit(
        x = z[train_idx, , drop = FALSE],
        y = keras::to_categorical(
            clusters[train_idx] - 1,
            num_classes = length(uclusters)
        ),
        epochs = epochs,
        batch_size = batch_size,
        validation_split = if (length(data_split) >= 2)
            data_split[2] / (data_split[1] + data_split[2]) else 0.2,
        verbose = 1
    )

    convt <- innsight::Converter$new(
        model,
        input_names = colnames(z),
        output_names = uclusters
    )
    method <- innsight::LRP$new(convt, z[test_idx, , drop = FALSE])
    method$get_result(type = "data.frame")
}

# ============================================================================
# Old code — kept as comment for reference
# ============================================================================

# #' Train weights for the markers using multinomial logistic regression
# #'
# #' @importFrom nnet multinom
# #' @importFrom tidyr pivot_longer
# #'
# #' @param path_to_gs Path to the gene set file without weights
# #' @param exprs The expression matrix, or a seurat object
# #'  (rows: genes, columns: samples/cells)
# #' @param level The level of the gene sets to train weights for if
# #'  you have multiple levels of gene sets.
# #' @param scaled Whether the expression matrix is scaled
# #' @param clusters A named vector of cluster ids
# #'  If `exprs` is a seurat object, this is ignored. The cluster ids are
# #'  taken from the seurat object.
# #' @param range The range of the weights
# #' @param data_split A vector of fractions for training and
# #'  testing. If only 1 fraction are provided, no testing set will be used.
# #' @param run_weights_on_test Whether to run the weights on the test set.
# #'  Requires that `data_split` has 2 elements.
# #'
# #' @return A data frame with the weights, that can be used directly by
# #'  \\code{\\link{gs_prepare}}.
# #'
# #' @export
# train_weights_mlr <- function(
#     path_to_gs,
#     exprs,
#     level = 1,
#     scaled = FALSE,
#     clusters = NULL,
#     range = c(1, 5),
#     data_split = c(0.7, 0.3),
#     run_weights_on_test = TRUE
# ) {
#     set.seed(1)
#     data <- prepare_data_for_training(
#         path_to_gs,
#         exprs,
#         level,
#         scaled,
#         clusters
#     )
#
#     z <- data$z
#     z$cluster <- clusters[rownames(z)]
#
#     test_data <- NULL
#     if (length(data_split) == 2) {
#         test_idx <- sample(
#             seq_len(nrow(data$z)), floor(nrow(data$z) * data_split[2])
#         )
#         test_data <- data$z[test_idx, , drop = FALSE]
#         rest_idx <- setdiff(seq_len(nrow(data$z)), test_idx)
#         z <- z[rest_idx, , drop = FALSE]
#         clusters <- clusters[rest_idx]
#     }
#
#     model <- multinom(cluster ~ ., data = z)
#     result <- summary(model)$coefficients
#     result <- as.data.frame(result[, -1, drop = FALSE])
#     result$output_node <- rownames(result)
#     result <- result %>%
#         pivot_longer(
#             -"output_node",
#             names_to = "feature",
#             values_to = "value"
#         )
#
#     weights <- compile_weights(result, data$gs, level, range)
#     if (!is.null(test_data) && run_weights_on_test) {
#         run_weights_on_test_data(
#             weights,
#             exprs,
#             clusters,
#             scaled,
#             rownames(test_data_x)
#         )
#     }
#     weights
# }

# ============================================================================
# Shared helpers
# ============================================================================

#' Run compiled weights on test data
#'
#' @keywords internal
#'
#' @param db The weights data frame
#' @param exprs The expression matrix, or a seurat object
#' (rows: genes, columns: samples/cells)
#' @param clusters A named vector of cluster ids
#' If `exprs` is a seurat object, this is ignored. The cluster ids are
#' taken from the seurat object.
#' @param scaled Whether the expression matrix is scaled
#' @param test_data_idx The row names of the test data
run_weights_on_test_data <- function(
    db,
    exprs,
    clusters,
    scaled,
    test_data_idx
) {
    if ("Seurat" %in% class(exprs)) {
        clusters <- Idents(exprs)
        exprs <- Seurat::GetAssayData(exprs, layer = "data")
        scaled <- FALSE
    }

    db$level <- NULL
    gs <- gs_prepare(db)
    scores <- hitype_score(
        exprs[, test_data_idx, drop = FALSE],
        gs,
        scaled = scaled,
        # Scoring with learned weights: marker sensitivity double-penalizes
        # shared markers (see hitype_score()). When every cell type shares
        # the same pooled marker set, all sensitivities are 0 and the whole
        # score matrix is zero, making this on-test assignment meaningless.
        use_sensitivity = FALSE
    )
    types <- hitype_assign(
        clusters[test_data_idx],
        scores[[1]],
        gs,
        threshold = 0
    )
    print(types)
}

#' Prepare the data for weight training
#'
#' @keywords internal
#'
#' @importFrom Seurat Idents
#'
#' @param path_to_gs Path to the gene set file without weights
#' @param exprs The expression matrix, or a seurat object
#'  (rows: genes, columns: samples/cells)
#' @param level The level of the gene sets to train weights for if
#'  you have multiple levels of gene sets.
#' @param scaled Whether the expression matrix is scaled
#' @param clusters The cell types of the cells. When `exprs` is a Seurat
#'  object, it can be `NULL` (default) to take the cell types from
#'  `Seurat::Idents()`, or a column name in the `meta.data` of the Seurat
#'  object that holds the cell type of each cell. When `exprs` is a matrix,
#'  it should be a named vector of cell types (names - cell names).
#'
#' @return A list with the gene sets, the z matrix and the clusters
prepare_data_for_training <- function(
    path_to_gs,
    exprs,
    level = 1,
    scaled = FALSE,
    clusters = NULL
) {
    if ("Seurat" %in% class(exprs)) {
        if (is.null(clusters)) {
            clusters <- Idents(exprs)
        } else if (is.character(clusters) && length(clusters) == 1) {
            if (!clusters %in% colnames(exprs@meta.data)) {
                stop(
                    paste0(
                        "The column `", clusters,
                        "` does not exist in the `meta.data` of the ",
                        "Seurat object. Give a column name that holds the ",
                        "cell type of each cell, or `NULL` to use ",
                        "`Seurat::Idents()`."
                    )
                )
            }
            # Named by the cells, like `Idents()`
            clusters <- exprs@meta.data[[clusters]]
            names(clusters) <- rownames(exprs@meta.data)
        } else {
            stop(
                paste0(
                    "When `exprs` is a Seurat object, `clusters` should be ",
                    "`NULL` (use `Seurat::Idents()`) or a single column ",
                    "name in the `meta.data` that holds the cell type of ",
                    "each cell."
                )
            )
        }
        exprs <- Seurat::GetAssayData(exprs, layer = "data")
        scaled <- FALSE
    }

    if (is.null(clusters)) {
        stop("Please provide a named vector of cluster ids")
    }

    gs <- gs_prepare(path_to_gs)$gene_sets[[level]]
    cell_types <- unique(as.character(clusters))
    file_types <- names(gs)
    has_cells <- file_types %in% cell_types

    # The training cell types come from the data: a cell type of the marker
    # file with no cells of that type in the data cannot be trained. It is
    # ignored (with a warning) and its markers are pooled for the cell types
    # of the data that are not covered by the marker file, so that those
    # cells are not left without markers.
    if (any(!has_cells)) {
        warning(
            paste0(
                "The following cell types in the marker file have no cells ",
                "of that type in the data and are ignored (their markers ",
                "are pooled for cell types not covered by the marker ",
                "file): ",
                paste(file_types[!has_cells], collapse = ", ")
            ),
            immediate. = TRUE
        )
    }
    leftover_markers <- unique(unlist(
        lapply(gs[!has_cells], function(x) x$markers)
    ))
    gs <- gs[has_cells]
    rest <- setdiff(cell_types, file_types)
    if (length(rest) > 0 && length(leftover_markers) > 0) {
        pooled <- rep(
            list(list(
                markers = leftover_markers,
                weights = rep(1, length(leftover_markers))
            )),
            length(rest)
        )
        names(pooled) <- rest
        gs <- c(gs, pooled)
    }
    if (length(gs) == 0) {
        stop("No markers are available for training.")
    }
    untrained <- setdiff(cell_types, names(gs))
    if (length(untrained) > 0) {
        warning(
            paste0(
                "The following cell types in the data are not covered by ",
                "the marker file and no markers are left for them; they ",
                "are not trained: ",
                paste(untrained, collapse = ", ")
            ),
            immediate. = TRUE
        )
    }

    all_markers <- na.omit(unlist(lapply(gs, function(x) x$markers)))
    non_exist_markers <- setdiff(all_markers, rownames(exprs))
    if (length(non_exist_markers) > 0) {
        warning(
            paste(
                "The following markers do not exist in the expression matrix:",
                paste(non_exist_markers, collapse = ", ")
            ),
            immediate. = TRUE
        )
    }
    all_markers <- setdiff(all_markers, non_exist_markers)
    for (ct in names(gs)) {
        gs[[ct]]$markers <- intersect(gs[[ct]]$markers, all_markers)
    }
    # All the cells are kept in the training data; a cell type with no
    # markers still contributes to the training of the other types as part
    # of the "rest" class.
    exprs <- exprs[all_markers, names(clusters), drop = FALSE]
    if (any(class(exprs) %in% c("dgCMatrix", "dgTMatrix"))) {
        exprs <- Matrix::t(exprs)
    } else {
        exprs <- t(exprs)
    }
    # The per-gene centering/scaling applied to the training matrix is
    # recorded so that model-based scoring (hitype_score_models()) can
    # reproduce the training values exactly. The statistics describe the
    # matrix before the test-set rows are removed in train_weights(): the
    # models are fit on rows scaled with these statistics, so scoring must
    # use the same ones — recomputing them on the training rows only would
    # not reproduce the values the models were fit on.
    if (!scaled) {
        z <- scale(exprs)
        center <- attr(z, "scaled:center")
        sds <- attr(z, "scaled:scale")
        # Genes with no variance across the cells (e.g. unexpressed) become
        # all-NA columns after scaling. They carry no signal for training,
        # yet make every model fit fail with "x has missing values" —
        # silently, when the fits are wrapped in tryCatch (all-zero
        # coefficients -> the rescale midpoint as weights). Drop them from
        # the matrix and the gene sets so they never reach the model or the
        # compiled weights.
        zero_var <- Matrix::colSums(is.na(z)) > 0
        if (any(zero_var)) {
            dropped <- colnames(z)[zero_var]
            warning(
                paste0(
                    "The following markers have no variance across the ",
                    "cells and are dropped: ",
                    paste(dropped, collapse = ", ")
                ),
                immediate. = TRUE
            )
            z <- z[, !zero_var, drop = FALSE]
            center <- center[colnames(z)]
            sds <- sds[colnames(z)]
            for (ct in names(gs)) {
                gs[[ct]]$markers <- intersect(
                    gs[[ct]]$markers, colnames(z)
                )
            }
        }
        exprs <- z
    } else {
        # Already-scaled input is used as-is: no centering/scaling was
        # applied to the training matrix, so scoring must not apply any
        # either (the bundle records scaled_input = TRUE for that)
        center <- setNames(rep(0, ncol(exprs)), colnames(exprs))
        sds <- setNames(rep(1, ncol(exprs)), colnames(exprs))
    }

    list(
        gs = gs, z = exprs, clusters = clusters,
        center = center, scale = sds, scaled_input = scaled
    )
}

#' Compile the weights
#'
#' @keywords internal
#'
#' @param weights A data frame with the weights
#' @param gs The gene sets
#' @param level The level of the gene sets
#' @param range The quantization range for `format = "db"` only: the
#'   4-element form `c(-low, -high, low, high)` (default `c(-5, -1, 1, 5)`)
#'   scales the positive weights into `[low, high]` and the negative ones
#'   into `[-high, -low]` before rounding them to integers (the db format
#'   can only encode integer weights via `+`/`-` repeats). The universal
#'   output is never rescaled — it keeps the raw trained weights.
#' @param format The format of the output data frame. One of
#'   `"universal"` (default) or `"db"`.
#' @param drop_zero Whether to drop markers whose trained weight is exactly
#'   zero from the output (default `TRUE`). Zero weights mean the marker
#'   carries no learned direction for the cell type; dropping them (before
#'   any db quantization) makes zero mean "no direction" — the marker is
#'   absent from the cell type's list. Pass `FALSE` to keep every candidate
#'   marker.
#' @param pos_only Whether to keep only the markers whose trained weight
#'   is positive in the output (default `FALSE`). A marker trained with a
#'   negative weight anti-correlates with the cell type and is reported with
#'   `direction = "negative"`; pass `TRUE` to drop those rows (as well as any
#'   exactly-zero ones, so `drop_zero` is then redundant) and keep only the
#'   markers overexpressed in the cell type.
#'
#' @return A data frame with the compiled weights in the universal marker
#'   format (default) or the db format (`format = "db"`), consumable by
#'   [gs_prepare()].
compile_weights <- function(
    weights, gs, level, range = c(-5, -1, 1, 5),
    format = c("universal", "db"), drop_zero = TRUE, pos_only = FALSE
) {
    format <- match.arg(format)
    if (format == "universal") {
        if (!identical(range, c(-5, -1, 1, 5))) {
            warning(
                "`range` only applies when `format = \"db\"`; the universal ",
                "output keeps the raw trained weights as-is."
            )
        }
    } else {
        if (
            length(range) != 4 || any(diff(range) < 0) ||
            range[2] > 0 || range[3] < 0
        ) {
            stop(
                "With `format = \"db\"`, `range` must have the 4-element ",
                "form c(-low, -high, low, high) (e.g. c(-5, -1, 1, 5)) so ",
                "that positive weights always rescale into the positive ",
                "band and negative ones into the negative band. A 2-element ",
                "range like c(-5, 5) would turn genuinely positive markers ",
                "negative."
            )
        }
    }
    weights <- weights %>%
        group_by(output_node, feature) %>%
        summarise(weight = mean(value), .groups = "drop")
    if (drop_zero) {
        # zero = no direction: drop exactly-zero weights (before any db
        # quantization) so a marker in the output is either positive or
        # negative
        weights <- weights[weights$weight != 0, , drop = FALSE]
    }
    if (pos_only) {
        # Keep only the markers whose learned direction is positive: a
        # negative weight (reported as `direction = "negative"`) or a zero
        # weight carries no "this marker marks the type" signal
        weights <- weights[weights$weight > 0, , drop = FALSE]
        # positive-only subsumes the zero drop, and markers filtered out
        # here must also drop out of each type's marker list below
        drop_zero <- TRUE
    }

    if (format == "universal") {
        # Raw weights: no rescaling, the trained coefficients are the
        # weights. gs_prepare() turns (direction, abs(weight)) back into
        # the signed coefficient, which is what hitype_score() consumes.
        parts <- lapply(names(gs), function(x) {
            markers <- explode(gs[[x]]$markers)
            w <- weights$weight[weights$output_node == x]
            names(w) <- weights$feature[weights$output_node == x]
            # name-aligned lookup: markers without a trained weight give NA
            v <- unname(w[markers])
            if (drop_zero) {
                # markers with zero (or no) trained weight have no direction:
                # drop them so absence in the list means "no direction"
                keep <- !is.na(v)
                markers <- markers[keep]
                v <- v[keep]
                if (length(markers) == 0) {
                    warning(
                        paste0(
                            "All markers of cell type '", x, "' have ",
                            if (pos_only) {
                                "no positive trained weight"
                            } else {
                                "zero or missing trained weights"
                            },
                            "; the cell type is dropped from the output"
                        ),
                        immediate. = TRUE
                    )
                }
            } else {
                v[is.na(v)] <- 1
            }
            data.frame(
                cell_type = x,
                gene = markers,
                direction = ifelse(v < 0, "negative", "positive"),
                weight = abs(v),
                level = as.integer(level),
                stringsAsFactors = FALSE
            )
        })
        return(do.call(rbind, parts))
    }

    # db format: quantize the raw weights into the signed integer bands of
    # `range`, then sign-encode them so they round-trip exactly through
    # gs_prepare()'s decoder (bare gene = 1, "+" x n = n + 1, "-" x n = -n,
    # "*" = 0). A slice with a single distinct value (or none) carries no
    # relative information; put it at the band midpoint.
    pos <- weights$weight > 0
    neg <- weights$weight < 0
    quantize_band <- function(v, to) {
        if (length(v) == 0 || length(unique(v)) < 2) {
            rep(round(mean(to)), length(v))
        } else {
            round(scales::rescale(v, to = to))
        }
    }
    weights$weight[pos] <- quantize_band(weights$weight[pos], range[3:4])
    weights$weight[neg] <- quantize_band(weights$weight[neg], range[1:2])

    db <- data.frame(
        cellName = names(gs),
        level = level,
        geneSymbolmore2 = rep("", length(gs))
    )
    db$geneSymbolmore1 <- unlist(lapply(names(gs), function(x) {
        markers <- explode(gs[[x]]$markers)
        w <- weights$weight[weights$output_node == x]
        names(w) <- weights$feature[weights$output_node == x]
        # name-aligned lookup, same semantics as the universal branch
        weight <- unname(w[markers])
        if (drop_zero) {
            keep <- !is.na(weight)
            markers <- markers[keep]
            weight <- weight[keep]
            if (length(markers) == 0) {
                warning(
                    paste0(
                        "All markers of cell type '", x, "' have zero or ",
                        "missing trained weights; the cell type is ",
                        "dropped from the output"
                    ),
                    immediate. = TRUE
                )
            }
        } else {
            weight[is.na(weight)] <- 1
        }
        markers <- sapply(seq_along(markers), function(i) {
            wi <- weight[i]
            if (wi > 0) {
                # "+" x (w - 1): the decoder adds 1 to the repeat count
                suffix <- paste0(rep("+", wi - 1), collapse = "")
            } else if (wi < 0) {
                suffix <- paste0(rep("-", -wi), collapse = "")
            } else {
                suffix <- "*"
            }
            paste0(markers[i], suffix, collapse = "")
        })
        paste0(markers, collapse = ",")
    }))
    db
}
