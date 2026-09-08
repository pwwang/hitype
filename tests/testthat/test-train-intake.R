# train_weights() intake: the training cell types come from the data
# (`clusters`), not from the marker file. A cell type in the marker file
# with no cells of that type in the data is ignored (with a warning) and
# its markers are pooled for the cell types in the data not covered by the
# marker file. No cells are dropped from the training data.

set.seed(123)
file_markers <- data.frame(
    cell_type = c("CT1", "CT1", "CT1", "CT2", "CT2"),
    gene = c("G1", "G2", "G3", "G4", "G5"),
    stringsAsFactors = FALSE
)
exprs <- matrix(
    runif(5 * 12, 0.1, 2),
    nrow = 5,
    dimnames = list(paste0("G", 1:5), paste0("c", 1:12))
)
# Named per-cell types from per-type cell counts, e.g. c(CT1 = 6, CT2 = 6)
types_of <- function(tab) {
    types <- rep(names(tab), tab)
    setNames(types, paste0("c", seq_along(types)))
}
aligned <- function(data) {
    # z rows (cells) and the returned clusters stay positionally aligned
    expect_equal(rownames(data$z), names(data$clusters))
    expect_equal(length(data$clusters), nrow(data$z))
}

test_that("fully matched marker file: own markers, silent, all cells kept", {
    clusters <- types_of(c(CT1 = 6, CT2 = 6))
    expect_silent(
        data <- hitype:::prepare_data_for_training(file_markers, exprs,
            clusters = clusters)
    )
    expect_equal(names(data$gs), c("CT1", "CT2"))
    expect_equal(data$gs$CT1$markers, c("G1", "G2", "G3"))
    expect_equal(data$gs$CT2$markers, c("G4", "G5"))
    expect_equal(colnames(data$z), c("G1", "G2", "G3", "G4", "G5"))
    aligned(data)
})

test_that("marker type without cells: ignored w/ warning, markers pooled", {
    # The user's example: file = CT1{G1..G3} + CT2{G4,G5}, data = CT1 + CT3
    clusters <- types_of(c(CT1 = 8, CT3 = 4))
    expect_warning(
        data <- hitype:::prepare_data_for_training(file_markers, exprs,
            clusters = clusters),
        "ignored.*CT2"
    )
    # CT1 keeps its own markers; CT3 is trained on the pool of CT2's
    # markers only (G4, G5) — not on CT1's markers
    expect_equal(names(data$gs), c("CT1", "CT3"))
    expect_equal(data$gs$CT1$markers, c("G1", "G2", "G3"))
    expect_equal(data$gs$CT3$markers, c("G4", "G5"))
    aligned(data)
})

test_that("end-to-end: pooled type gets output rows, ignored type not", {
    clusters <- types_of(c(CT1 = 8, CT3 = 4))
    expect_warning(
        w <- train_weights(
            path_to_gs = file_markers,
            exprs = exprs,
            clusters = clusters,
            method = "uniform",
            data_split = c(1)
        ),
        "ignored.*CT2"
    )
    # No CT2 is reported; CT1 rows over its own markers, CT3 rows over the
    # pooled CT2 markers
    expect_equal(sort(unique(w$cell_type)), c("CT1", "CT3"))
    expect_setequal(w$gene[w$cell_type == "CT1"], c("G1", "G2", "G3"))
    expect_setequal(w$gene[w$cell_type == "CT3"], c("G4", "G5"))
})

test_that("data type without file entry nor pooled markers is not trained", {
    clusters <- types_of(c(CT1 = 4, CT2 = 4, CT3 = 4))
    expect_warning(
        data <- hitype:::prepare_data_for_training(file_markers, exprs,
            clusters = clusters),
        "not covered.*CT3"
    )
    expect_equal(names(data$gs), c("CT1", "CT2"))
    # The unmatched cells stay in the training data ("rest" class)
    aligned(data)
})

test_that("zero overlap: all data types train on the pooled markers", {
    clusters <- types_of(c(CT3 = 6, CT4 = 6))
    expect_warning(
        data <- hitype:::prepare_data_for_training(file_markers, exprs,
            clusters = clusters),
        "ignored.*CT1, CT2"
    )
    expect_equal(names(data$gs), c("CT3", "CT4"))
    expect_equal(data$gs$CT3$markers, c("G1", "G2", "G3", "G4", "G5"))
    expect_equal(data$gs$CT4$markers, c("G1", "G2", "G3", "G4", "G5"))
    aligned(data)
})

test_that("Seurat input: clusters is NULL (Idents) or a meta.data column", {
    skip_if_not_installed("Seurat")
    obj <- Seurat::CreateSeuratObject(
        counts = Matrix::Matrix(round(exprs * 10), sparse = TRUE),
        min.cells = 0, min.features = 0
    )
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    types <- factor(rep(c("CT1", "CT2"), each = 6), levels = c("CT1", "CT2"))
    names(types) <- paste0("c", 1:12)
    Seurat::Idents(obj) <- types
    obj@meta.data$cell_type <- as.character(types)

    # NULL -> the current Idents() provide the cell types
    expect_silent(
        data_null <- hitype:::prepare_data_for_training(file_markers, obj)
    )
    expect_equal(as.character(data_null$clusters), as.character(types))
    aligned(data_null)

    # A meta.data column name -> that column provides the cell types
    expect_silent(
        data_col <- hitype:::prepare_data_for_training(file_markers, obj,
            clusters = "cell_type")
    )
    expect_equal(as.character(data_col$clusters), as.character(types))
    aligned(data_col)

    # A non-existent column or a real clusters vector -> errors
    expect_error(
        hitype:::prepare_data_for_training(file_markers, obj,
            clusters = "not_a_column"),
        "does not exist in the `meta.data`"
    )
    expect_error(
        hitype:::prepare_data_for_training(file_markers, obj,
            clusters = c("CT1", "CT2")),
        "column name"
    )
})
