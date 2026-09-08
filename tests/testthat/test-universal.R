# Universal marker format (long table) intake and output
#
# A table is in the universal format when it has (case-insensitively,
# aliases allowed) a cell_type-family and a gene-family column.

# Read the decoded weights of a cell type as a named vector (names = genes)
gs_weights <- function(gs, level = 1, ct) {
    setNames(gs$gene_sets[[level]][[ct]]$weights, gs$gene_sets[[level]][[ct]]$markers)
}

test_that("universal tables are detected with aliases, case-insensitively", {
    markers <- data.frame(
        CellType = c("CD4 T", "CD4 T", "CD8 T", "CD8 T"),
        Marker = c("CD3E", "IL7R", "CD8A", "GZMB"),
        sign = c("pos", "Negative", "+", "-"),
        weight = c(2, 0.5, 3, 2.5)
    )
    gs <- gs_prepare(markers)
    expect_null(gs$cell_names)
    w <- gs_weights(gs, 1, "CD4 T")
    expect_equal(w, c(CD3E = 2, IL7R = -0.5))
    w <- gs_weights(gs, 1, "CD8 T")
    expect_equal(w, c(CD8A = 3, GZMB = -2.5))
})

test_that("universal decode rule: direction is authoritative", {
    # Direction positive/negative -> abs(weight) / -abs(weight), even when
    # the signed weight says otherwise; no direction but a weight column ->
    # the signed weight as-is; neither -> 1.
    markers <- data.frame(
        cell_type = c("A", "A", "A", "A", "A", "A", "A", "A"),
        gene = letters[1:8],
        direction = c(
            "positive", "negative", "positive", "negative",
            NA, "", NA, NA
        ),
        weight = c(1, 1, -3, -2.37, -3, 0.5, NA, NA),
        stringsAsFactors = FALSE
    )
    gs <- gs_prepare(markers)
    expect_equal(
        gs_weights(gs, 1, "A"),
        c(a = 1, b = -1, c = 3, d = -2.37, e = -3, f = 0.5, g = 1, h = 1)
    )
})

test_that("universal decode keeps fractional weights exact", {
    markers <- data.frame(
        cell_type = c("A", "A"),
        gene = c("x", "y"),
        direction = c("positive", "negative"),
        weight = c(1.23456789, 0.000123)
    )
    gs <- gs_prepare(markers)
    expect_equal(
        gs_weights(gs, 1, "A"),
        c(x = 1.23456789, y = -0.000123)
    )
})

test_that("universal direction aliases are accepted, invalid values error", {
    ok <- data.frame(
        cell_type = "A", gene = "x",
        direction = c("positive", "pos", "+", "POS", "Positive",
                      "negative", "neg", "-", "NEG")
    )
    expect_silent(gs_prepare(ok))
    bad <- data.frame(cell_type = "A", gene = "x", direction = "up")
    expect_error(gs_prepare(bad), "Accepted values")
    bad2 <- data.frame(cell_type = "A", gene = "x", direction = "pos;up")
    expect_error(gs_prepare(bad2), "pos;up")
})

test_that("universal non-numeric weights error", {
    bad <- data.frame(
        cell_type = "A", gene = "x", direction = "positive", weight = "heavy"
    )
    expect_error(gs_prepare(bad), "non-numeric")
})

test_that("universal rows with NA/empty cell types or genes are dropped", {
    markers <- data.frame(
        cell_type = c("A", "A", NA, " ", "A"),
        gene = c("x", "y", "z", "w", NA)
    )
    gs <- gs_prepare(markers)
    expect_equal(names(gs$gene_sets[[1]]), "A")
    expect_equal(gs$gene_sets[[1]]$A$markers, c("x", "y"))
})

test_that("universal duplicate genes within a cell type keep the first", {
    markers <- data.frame(
        cell_type = c("A", "A", "A"),
        gene = c("x", "x", "y"),
        weight = c(2, 3, 4)
    )
    gs <- gs_prepare(markers)
    expect_equal(
        gs_weights(gs, 1, "A"),
        c(x = 2, y = 4)
    )
})

test_that("universal cell types must not contain `,` or `;`", {
    markers <- data.frame(cell_type = "A,B", gene = "x")
    expect_error(gs_prepare(markers), "The cell names should not contain")
})

test_that("universal level defaults to 1 and is validated", {
    m1 <- data.frame(cell_type = "A", gene = "x")
    expect_length(gs_prepare(m1)$gene_sets, 1)

    m2 <- data.frame(cell_type = "A", gene = "x", level = 2)
    expect_error(gs_prepare(m2), "Level should start from 1.")

    m3 <- data.frame(
        cell_type = c("A", "A"), gene = c("x", "y"), level = c(1, 3)
    )
    expect_error(gs_prepare(m3), "Level should be consecutive.")

    # Multi-level: one gene_sets entry per level, no cell_names hierarchy
    m4 <- data.frame(
        cell_type = c("A", "A", "B"), gene = c("x", "y", "z"),
        level = c(1, 1, 2)
    )
    gs <- gs_prepare(m4)
    expect_null(gs$cell_names)
    expect_equal(names(gs$gene_sets), c("1", "2"))
    expect_equal(gs$gene_sets[["1"]]$A$markers, c("x", "y"))
    expect_equal(gs$gene_sets[["2"]]$B$markers, "z")
})

test_that("universal tissue filter works with `tissue` and its alias", {
    markers <- data.frame(
        cell_type = c("A", "B"),
        gene = c("x", "y"),
        tissue = c("blood", "liver")
    )
    gs <- gs_prepare(markers, tissue_type = "liver")
    expect_equal(names(gs$gene_sets[[1]]), "B")

    markers2 <- data.frame(
        cell_type = c("A", "B"),
        gene = c("x", "y"),
        tissueType = c("blood", "liver")
    )
    gs2 <- gs_prepare(markers2, tissue_type = "blood")
    expect_equal(names(gs2$gene_sets[[1]]), "A")

    no_tissue <- data.frame(cell_type = c("A", "B"), gene = c("x", "y"))
    expect_error(
        gs_prepare(no_tissue, tissue_type = "blood"),
        "does not have the `tissue` column"
    )
})

test_that("universal weight_encoding is applied", {
    markers <- data.frame(
        cell_type = c("A", "A"),
        gene = c("x", "y"),
        direction = c("positive", "negative"),
        weight = c(2, 1)
    )
    gs <- gs_prepare(markers, weight_encoding = function(x) x * x)
    expect_equal(
        gs_weights(gs, 1, "A"),
        c(x = 4, y = 1)
    )
})

test_that("universal marker files (.csv/.rds) are auto-detected", {
    markers <- data.frame(
        cell_type = c("A", "A", "B"),
        gene = c("x", "y", "z"),
        direction = c("positive", "negative", "positive"),
        weight = c(2, 1, 3)
    )
    csv <- tempfile(fileext = ".csv")
    write.csv(markers, csv, row.names = FALSE)
    gs <- gs_prepare(csv)
    expect_equal(
        gs_weights(gs, 1, "A"),
        c(x = 2, y = -1)
    )
    expect_equal(gs_weights(gs, 1, "B"), c(z = 3))

    rds <- tempfile(fileext = ".rds")
    saveRDS(markers, rds)
    gs2 <- gs_prepare(rds)
    expect_equal(
        gs_weights(gs2, 1, "A"),
        c(x = 2, y = -1)
    )
})

test_that("compile_weights() universal output round-trips exact weights", {
    gs <- list(
        `CD4 T` = list(markers = c("IL7R", "CCR7", "S100A4")),
        `CD8 T` = list(markers = c("CD8A", "CD8B", "GZMB"))
    )
    weights <- data.frame(
        output_node = c(
            "CD4 T", "CD4 T", "CD4 T",
            "CD8 T", "CD8 T", "CD8 T"
        ),
        feature = c(
            "IL7R", "CCR7", "S100A4",
            "CD8A", "CD8B", "GZMB"
        ),
        # Span exactly the range so the rescale below is the identity
        value = c(2, -1, 0, -2, 1, 1.5)
    )
    compiled <- compile_weights(weights, gs, level = 1, range = c(-2, 2))
    expect_equal(
        colnames(compiled),
        c("cell_type", "gene", "direction", "weight", "level")
    )
    expect_true(all(compiled$weight >= 0))
    expect_equal(
        compiled$direction[compiled$weight == 0],
        rep("positive", sum(compiled$weight == 0))
    )
    # Re-decoded through gs_prepare: exact signed weights, incl. zero
    gs2 <- gs_prepare(compiled)
    expect_null(gs2$cell_names)
    expect_equal(
        gs_weights(gs2, 1, "CD4 T"),
        c(IL7R = 2, CCR7 = -1, S100A4 = 0)
    )
    expect_equal(
        gs_weights(gs2, 1, "CD8 T"),
        c(CD8A = -2, CD8B = 1, GZMB = 1.5)
    )

    # Legacy db format is still available
    db <- compile_weights(weights, gs, level = 1, range = c(-2, 2),
        format = "db")
    expect_true(all(
        c("cellName", "geneSymbolmore1", "geneSymbolmore2", "level") %in%
            colnames(db)
    ))
    expect_equal(nrow(db), 2)
})

test_that("universal marker files (.qs/.qs2) are auto-detected", {
    markers <- data.frame(
        cell_type = c("A", "B"), gene = c("x", "y"),
        direction = "negative", weight = c(2, 1)
    )
    if (requireNamespace("qs2", quietly = TRUE)) {
        qfile <- tempfile(fileext = ".qs2")
        qs2::qs_save(markers, qfile)
        expect_equal(
            gs_weights(gs_prepare(qfile), 1, "A"),
            c(x = -2)
        )
    } else {
        expect_error(gs_prepare("x.qs2"), "package `qs2` is required")
    }
})
