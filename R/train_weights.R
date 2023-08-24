#' Train weights for the markers
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @import keras
#' @importFrom innsight Converter
#' @importFrom innsight LRP
#'
#' @param path_to_gs Path to the gene set file without weights
#' @param exprs The expression matrix, or a seurat object
#'  (rows: genes, columns: samples/cells)
#' @param level The level of the gene sets to train weights for if
#'  you have multiple levels of gene sets.
#' @param scaled Whether the expression matrix is scaled
#' @param clusters A named vector of cluster ids
#'  If `exprs` is a seurat object, this is ignored. The cluster ids are
#'  taken from the seurat object.
#' @param range The range of the weights
#' @param data_split A vector of fractions for training, validation and
#'  testing. If only two fractions are provided, no testing set will be used.
#' @param epochs The number of epochs to train
#' @param batch_size The batch size
#' @param run_weights_on_test Whether to run the weights on the test set.
#'  Requires that `data_split` has three elements.
#'
#' @return A data frame with the weights, that can be used directly by
#'  \code{\link{gs_prepare}}.
#'
#' @export
train_weights_nn <- function(
    path_to_gs,
    exprs,
    level = 1,
    scaled = FALSE,
    clusters = NULL,
    range = c(1, 5),
    data_split = c(0.7, 0.2, 0.1),
    epochs = 20,
    batch_size = 32,
    run_weights_on_test = TRUE
) {
    set.seed(1)
    data <- prepare_data_for_training(
        path_to_gs,
        exprs,
        level,
        scaled,
        clusters
    )

    clusters <- as.integer(data$clusters)
    uclusters <- unique(clusters)

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
    model <- keras_model_sequential() %>%
        layer_masking(mask_value = 0, input_shape = ncol(data$z)) %>%
        layer_dense(
            units = 64, activation = "relu", input_shape = ncol(data$z)
        ) %>%
        layer_dropout(rate = 0.2) %>%
        layer_dense(units = 64, activation = "relu") %>%
        layer_dropout(rate = 0.2) %>%
        layer_dense(units = length(uclusters), activation = "softmax")

    model %>% compile(
        loss = "categorical_crossentropy",
        optimizer = "adam",
        metrics = c("accuracy")
    )
    model %>% fit(
        x = data$z,
        y = to_categorical(clusters - 1, num_classes = length(uclusters)),
        epochs = epochs,
        batch_size = batch_size,
        validation_split = data_split[2] / (data_split[1] + data_split[2]),
        verbose = 1
    )

    if (!is.null(test_data_x)) {
        cat("Evaluating on test data\n")
        test_data_y <- to_categorical(
            test_data_y - 1,
            num_classes = length(uclusters)
        )
        model %>% evaluate(
            x = test_data_x,
            y = test_data_y,
            batch_size = batch_size,
            verbose = 1
        )
    }

    convt <- Converter$new(
        model,
        input_names = colnames(data$z),
        output_names = unique(data$clusters)
    )
    method <- LRP$new(convt, data$z)
    result <- method$get_result(type = "data.frame")

    weights <- compile_weights(result, data$gs, level, range)
    if (!is.null(test_data_x) && run_weights_on_test) {
        run_weights_on_test_data(
            weights,
            exprs,
            clusters,
            scaled,
            rownames(test_data_x)
        )
    }
    weights
}

#' Train weights for the markers using multinomial logistic regression
#'
#' @importFrom nnet multinom
#'
#' @param path_to_gs Path to the gene set file without weights
#' @param exprs The expression matrix, or a seurat object
#'  (rows: genes, columns: samples/cells)
#' @param level The level of the gene sets to train weights for if
#'  you have multiple levels of gene sets.
#' @param scaled Whether the expression matrix is scaled
#' @param clusters A named vector of cluster ids
#'  If `exprs` is a seurat object, this is ignored. The cluster ids are
#'  taken from the seurat object.
#' @param range The range of the weights
#' @param data_split A vector of fractions for training and
#'  testing. If only 1 fraction are provided, no testing set will be used.
#' @param run_weights_on_test Whether to run the weights on the test set.
#'  Requires that `data_split` has 2 elements.
#'
#' @return A data frame with the weights, that can be used directly by
#'  \code{\link{gs_prepare}}.
#'
#' @export
train_weights_mlr <- function(
    path_to_gs,
    exprs,
    level = 1,
    scaled = FALSE,
    clusters = NULL,
    range = c(1, 5),
    data_split = c(0.7, 0.3),
    run_weights_on_test = TRUE
) {
    set.seed(1)
    data <- prepare_data_for_training(
        path_to_gs,
        exprs,
        level,
        scaled,
        clusters
    )

    z <- data$z
    z$cluster <- clusters[rownames(z)]

    test_data <- NULL
    if (length(data_split) == 2) {
        test_idx <- sample(
            seq_len(nrow(data$z)), floor(nrow(data$z) * data_split[2])
        )
        test_data <- data$z[test_idx, , drop = FALSE]
        rest_idx <- setdiff(seq_len(nrow(data$z)), test_idx)
        z <- z[rest_idx, , drop = FALSE]
        clusters <- clusters[rest_idx]
    }

    model <- multinom(cluster ~ ., data = z)
    result <- summary(model)$coefficients
    result <- as.data.frame(result[, -1, drop = FALSE])
    result$output_node <- rownames(result)
    result <- result %>%
        tidyr::pivot_longer(
            -"output_node",
            names_to = "feature",
            values_to = "value"
        )

    weights <- compile_weights(result, data$gs, level, range)
    if (!is.null(test_data) && run_weights_on_test) {
        run_weights_on_test_data(
            weights,
            exprs,
            clusters,
            scaled,
            rownames(test_data_x)
        )
    }
    weights
}

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
        exprs <- exprs@assays$RNA@data
        scaled <- FALSE
    }

    db$level <- NULL
    gs <- gs_prepare(db)
    scores <- hitype_score(
        exprs[, test_data_idx, drop = FALSE], gs, scaled = scaled
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
#' @param clusters A named vector of cluster ids
#'  If `exprs` is a seurat object, this is ignored. The cluster ids are
#'  taken from the seurat object.
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
        clusters <- Idents(exprs)
        exprs <- exprs@assays$RNA@data
        scaled <- FALSE
    }

    if (is.null(clusters)) {
        stop("Please provide a named vector of cluster ids")
    }

    gs <- gs_prepare(path_to_gs)$gene_sets[[level]]
    non_exist_clusters <- setdiff(names(gs), unique(clusters))
    if (length(non_exist_clusters) > 0) {
        stop(
            paste(
                "The following clusters do not exist in the expression matrix:",
                paste(non_exist_clusters, collapse = ", ")
            )
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
    exprs <- exprs[
        all_markers,
        names(clusters[clusters %in% names(gs)]),
        drop = FALSE
    ]
    if (any(class(exprs) %in% c("dgCMatrix", "dgTMatrix"))) {
        exprs <- Matrix::t(exprs)
    } else {
        exprs <- t(exprs)
    }
    if (!scaled) { exprs <- scale(exprs) }

    mask <- matrix(FALSE, nrow = nrow(exprs), ncol = ncol(exprs))
    rownames(mask) <- rownames(exprs)
    colnames(mask) <- colnames(exprs)
    for (ct in names(gs)) {
        mask[
            rownames(mask) %in% names(clusters[clusters == ct]),
            gs[[ct]]$markers
        ] <- TRUE
    }
    exprs[!mask] <- 0
    # z : rows: samples, columns: genes
    list(gs = gs, z = exprs, clusters = clusters)
}

#' Compile the weights
#'
#' @keywords internal
#'
#' @param weights A data frame with the weights
#' @param gs The gene sets
#' @param level The level of the gene sets
#' @param range The range of the weights
#'
#' @return A data frame as the db
compile_weights <- function(weights, gs, level, range) {
    if (any(diff(range) < 0)) {
        stop("range must be increasing")
    }
    if (!length(range) %in% c(2, 4)) {
        stop("range must be of length 2 or 4")
    }
    if (length(range) == 4) {
        if (range[2] > 0 || range[3] < 0) {
            stop("range must be of the form c(-low, -high, low, high)")
        }
    }
    weights <- weights %>%
        group_by(output_node, feature) %>%
        summarise(weight = mean(value))
    if (length(range) == 2) {
        weights <- weights %>% mutate(weight = scales::rescale(weight, to = range))
    } else {
        weights$weight[weights$weight > 0] <- scales::rescale(
            weights$weight[weights$weight > 0],
            to = range[3:4]
        )
        weights$weight[weights$weight < 0] <- scales::rescale(
            weights$weight[weights$weight < 0],
            to = range[1:2]
        )
    }

    db <- data.frame(
        cellName = names(gs),
        level = level,
        geneSymbolmore2 = rep("", length(gs))
    )
    db$geneSymbolmore1 <- unlist(lapply(names(gs), function(x) {
        markers <- explode(gs[[x]]$markers)
        weight <- weights[
            weights$output_node == x & weights$feature %in% markers,
            "weight",
            drop = TRUE
        ]
        markers <- sapply(seq_along(markers), function(i) {
            if (weight[i] > 0) {
                sign <- "+"
            } else if (weight[i] < 0) {
                sign <- "-"
            } else {
                sign <- "*"
                weight[i] <- 1
            }
            suffix <- paste0(rep(sign, abs(weight[i])), collapse = "")
            paste0(markers[i], suffix, collapse = "")
        })
        paste0(markers, collapse = ",")
    }))
    db
}