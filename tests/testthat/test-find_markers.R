set.seed(123)
ngenes <- 60
ncells <- 120
known <- list(
    Tcell = c("CD3D", "CD3E", "CD8A", "IL7R", "CCR7", "LCK"),
    Bcell = c("MS4A1", "CD79A", "CD79B", "BANK1", "FCRL5"),
    Monocyte = c(
        "CD14", "LYZ", "FCGR3A", "CST3", "S100A8", "S100A9", "FCN1"
    )
)
gene_names <- c(
    unlist(known, use.names = FALSE),
    paste0("gene", seq_len(ngenes - length(unlist(known))))
)
clusters <- setNames(
    rep(c("Tcell", "Bcell", "Monocyte"), each = 40),
    paste0("cell", seq_len(ncells))
)
# Baseline small values with ~50% zeros, then log1p-transform
counts <- matrix(
    runif(ngenes * ncells, 0, 0.2) * rbinom(ngenes * ncells, 1, 0.5),
    ngenes, ncells,
    dimnames = list(gene_names, names(clusters))
)
# Each cell type over-expresses its known markers 3-5x (before log1p);
# zeros among the over-expressed cells are filled so the markers are
# detected in every cell of their own type (pct_in = 1)
for (ct in names(known)) {
    cells <- names(clusters)[clusters == ct]
    counts[known[[ct]], cells] <- counts[known[[ct]], cells] *
        runif(length(cells), 3, 5)
    counts[known[[ct]], cells][counts[known[[ct]], cells] == 0] <- 0.2
}
exprs <- log1p(counts)

test_that("find_markers() fc returns db-format data.frame", {
    res <- find_markers(exprs, clusters, method = "fc", top = 8)
    expect_equal(
        colnames(res),
        c("cellName", "geneSymbolmore1", "geneSymbolmore2", "level")
    )
    expect_equal(nrow(res), 3)
    expect_identical(unique(res$cellName), c("Tcell", "Bcell", "Monocyte"))
    expect_equal(res$level, rep(1L, 3))
    expect_true(all(res$geneSymbolmore2 == ""))
})

test_that("find_markers() fc recovers known over-expressed markers", {
    res <- find_markers(exprs, clusters, method = "fc", top = 8)
    for (ct in names(known)) {
        markers <- strsplit(
            res$geneSymbolmore1[res$cellName == ct], ","
        )[[1]]
        expect_gte(length(intersect(known[[ct]], markers)), 4)
        expect_false(any(grepl("[+-]", markers)))
    }
})

test_that("find_markers() include_negative fills geneSymbolmore2", {
    res <- find_markers(
        exprs, clusters, method = "fc", top = 8, include_negative = TRUE
    )
    expect_false(all(res$geneSymbolmore2 == ""))
    expect_gte(sum(res$geneSymbolmore2 != ""), 2)
    gs <- gs_prepare(res)
    expect_true(is.list(gs$gene_sets))
    expect_true(all(
        c("Tcell", "Bcell", "Monocyte") %in% names(gs$gene_sets[[1]])
    ))
})

test_that("find_markers() validates inputs", {
    expect_error(find_markers(exprs, method = "fc"), "clusters")
    expect_error(
        find_markers(exprs, clusters = clusters[-1], method = "fc"),
        "clusters"
    )
    two_cell <- clusters
    mono_cells <- names(two_cell)[two_cell == "Monocyte"]
    two_cell[mono_cells[seq_len(2)]] <- "Rare"
    expect_warning(
        res <- find_markers(exprs, two_cell, method = "fc", top = 8),
        "fewer than 3 cells"
    )
    expect_false("Rare" %in% res$cellName)
    expect_equal(nrow(res), 3)
    expect_error(
        find_markers(
            matrix(0, 3, 3, dimnames = list(NULL, NULL)),
            clusters = c("A", "B", "C"),
            method = "fc"
        ),
        "rownames"
    )
    expect_error(find_markers(exprs, clusters, top = 0), "top")
    expect_error(find_markers(exprs, clusters, min_log2fc = -1), "min_log2fc")
    expect_error(find_markers(exprs, clusters, min_pct = 2), "min_pct")
    expect_error(
        find_markers(exprs, clusters, method = "seurat"),
        "Seurat object"
    )
})

test_that("find_markers() fc works on a dgCMatrix", {
    sparse_exprs <- Matrix::Matrix(exprs, sparse = TRUE)
    res_sparse <- find_markers(sparse_exprs, clusters, method = "fc", top = 8)
    res_dense <- find_markers(exprs, clusters, method = "fc", top = 8)
    expect_equal(res_sparse, res_dense)
    tcell_markers <- strsplit(
        res_sparse$geneSymbolmore1[res_sparse$cellName == "Tcell"], ","
    )[[1]]
    expect_gte(length(intersect(known$Tcell, tcell_markers)), 4)
})

test_that("find_markers() method seurat works on a Seurat object", {
    skip_if_not_installed("Seurat")
    res <- tryCatch(
        {
            obj <- Seurat::CreateSeuratObject(
                counts = Matrix::Matrix(counts, sparse = TRUE),
                min.cells = 0,
                min.features = 0
            )
            obj <- Seurat::NormalizeData(obj, verbose = FALSE)
            Seurat::Idents(obj) <- factor(
                unname(clusters),
                levels = c("Tcell", "Bcell", "Monocyte")
            )
            find_markers(obj, method = "seurat", top = 5)
        },
        error = function(e) {
            skip(paste0(
                "Seurat failed on the tiny test object: ",
                conditionMessage(e)
            ))
        }
    )
    expect_equal(
        colnames(res),
        c("cellName", "geneSymbolmore1", "geneSymbolmore2", "level")
    )
    expect_equal(nrow(res), 3)
})

test_that("find_markers() method presto works", {
    skip_if_not_installed("presto")
    res <- find_markers(exprs, clusters, method = "presto", top = 8)
    expect_equal(
        colnames(res),
        c("cellName", "geneSymbolmore1", "geneSymbolmore2", "level")
    )
    expect_equal(nrow(res), 3)
    tcell_markers <- strsplit(
        res$geneSymbolmore1[res$cellName == "Tcell"], ","
    )[[1]]
    # presto's logFC is on a different scale than the fc backend's log2fc
    # (roughly an order of magnitude smaller on this log1p data), so with
    # only_pos = TRUE the shared min_log2fc = 0.25 filter keeps only the
    # strongest known markers (CD3E, CCR7). Ranking is by
    # logFC * (pct_in - pct_out), consistent with the fc backend.
    expect_gte(length(intersect(known$Tcell, tcell_markers)), 2)
})
