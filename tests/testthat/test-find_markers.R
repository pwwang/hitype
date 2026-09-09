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
    res <- find_markers(exprs, clusters, method = "fc", top = 8, format = "db")
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
    res <- find_markers(exprs, clusters, method = "fc", top = 8, format = "db")
    for (ct in names(known)) {
        markers <- strsplit(
            res$geneSymbolmore1[res$cellName == ct], ","
        )[[1]]
        expect_gte(length(intersect(known[[ct]], markers)), 4)
        expect_false(any(grepl("[+-]", markers)))
    }
})

test_that("find_markers() !pos_only fills geneSymbolmore2", {
    res <- find_markers(
        exprs, clusters, method = "fc", top = 8, pos_only = FALSE,
        format = "db"
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
        res <- find_markers(exprs, two_cell, method = "fc", top = 8,
            format = "db"),
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
    res_sparse <- find_markers(
        sparse_exprs, clusters, method = "fc", top = 8, format = "db"
    )
    res_dense <- find_markers(
        exprs, clusters, method = "fc", top = 8, format = "db"
    )
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
            find_markers(obj, method = "seurat", top = 5, format = "db")
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
    res <- find_markers(exprs, clusters, method = "presto", top = 8,
        format = "db")
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
    # pos_only = TRUE the shared min_log2fc = 0.25 filter keeps only the
    # strongest known markers (CD3E, CCR7). Ranking is by
    # logFC * (pct_in - pct_out), consistent with the fc backend.
    expect_gte(length(intersect(known$Tcell, tcell_markers)), 2)
})

test_that("find_markers() defaults to the universal marker format", {
    res <- find_markers(exprs, clusters, method = "fc", top = 8)
    expect_equal(
        colnames(res),
        c("cell_type", "gene", "direction", "level")
    )
    expect_gt(nrow(res), 3)  # one row per cell_type-gene pair
    expect_true(all(res$direction == "positive"))
    expect_equal(res$level, rep(1L, nrow(res)))
    expect_true(all(!is.na(res$cell_type) & !is.na(res$gene)))
    # Positive markers of each cell type are preserved
    for (ct in names(known)) {
        genes <- res$gene[res$cell_type == ct]
        expect_gte(length(intersect(known[[ct]], genes)), 4)
    }
})

test_that("find_markers() universal format carries negative markers", {
    res <- find_markers(
        exprs, clusters, method = "fc", top = 8, pos_only = FALSE
    )
    expect_true(any(res$direction == "negative"))
    expect_true(any(res$direction == "positive"))
    # Every cell type keeps its positive markers
    for (ct in names(known)) {
        genes <- res$gene[res$cell_type == ct & res$direction == "positive"]
        expect_gte(length(genes), 1)
    }
})

test_that("find_markers() universal output round-trips through gs_prepare", {
    res <- find_markers(
        exprs, clusters, method = "fc", top = 8, pos_only = FALSE
    )
    gs <- gs_prepare(res)
    expect_null(gs$cell_names)
    for (ct in c("Tcell", "Bcell", "Monocyte")) {
        markers <- gs$gene_sets[[1]][[ct]]$markers
        weights <- gs$gene_sets[[1]][[ct]]$weights
        names(weights) <- markers
        pos <- res$gene[res$cell_type == ct & res$direction == "positive"]
        neg <- res$gene[res$cell_type == ct & res$direction == "negative"]
        expect_true(all(weights[pos] == 1))
        expect_true(all(weights[neg] == -1))
    }
})

test_that("find_markers() validates format", {
    expect_error(
        find_markers(exprs, clusters, method = "fc", format = "wide"),
        "arg"
    )
})

# --- synthetic A/B/C fixtures for the `against` / max_pct_out tests ---
# A and B are correlated siblings sharing the pan-lineage genes (pan,
# pan2) that C lacks; disc separates A from B/C. Types are 20/20/40 cells.
sib_clusters <- setNames(
    rep(c("A", "B", "C"), c(20, 20, 40)),
    paste0("sibcell", seq_len(80))
)
sib_counts <- matrix(
    0, 6, 80,
    dimnames = list(
        c("pan", "pan2", "disc", "bg1", "bg2", "noise"),
        names(sib_clusters)
    )
)
sib_counts["pan", sib_clusters %in% c("A", "B")] <- 10
sib_counts["pan2", sib_clusters %in% c("A", "B")] <- 5
sib_counts["disc", sib_clusters == "A"] <- 10
sib_counts["disc", sib_clusters == "C"] <- 3
sib_counts["bg1", ] <- 1
sib_counts["bg2", ] <- 0.5
sib_counts["noise", sib_clusters == "A"][1] <- 0.5
exprs_sib <- log1p(sib_counts)

# collect the warnings raised by an expression (muffling them)
catch_warnings <- function(code) {
    out <- character(0)
    withCallingHandlers(
        code,
        warning = function(cond) {
            out <<- c(out, conditionMessage(cond))
            invokeRestart("muffleWarning")
        }
    )
    out
}

test_that("find_markers() fc `against` recovers the sibling discriminator", {
    res_null <- find_markers(
        exprs_sib, sib_clusters, method = "fc", top = 5, format = "db"
    )
    res_ag <- find_markers(
        exprs_sib, sib_clusters, method = "fc", top = 5,
        against = "B", format = "db"
    )
    a_null <- strsplit(
        res_null$geneSymbolmore1[res_null$cellName == "A"], ","
    )[[1]]
    a_ag <- strsplit(
        res_ag$geneSymbolmore1[res_ag$cellName == "A"], ","
    )[[1]]
    # the pan-lineage gene outranks the A-vs-B discriminator one-vs-rest
    expect_true(all(c("pan", "disc") %in% a_null))
    expect_lt(which(a_null == "pan")[1], which(a_null == "disc")[1])
    # ... but comparing A against its sibling B keeps only the discriminator
    expect_identical(a_ag, "disc")
    expect_false("pan" %in% a_ag)
    # the type listed in `against` is not compared against itself
    expect_identical(res_ag$geneSymbolmore1[res_ag$cellName == "B"], "")
})

test_that("find_markers() fc `against = 'nearest'` uses the most-correlated type", {
    # A's most-correlated type is B (mean-expression-profile Pearson
    # correlation: r(A,B) ~ 0.54 > r(A,C) ~ 0.27 on this fixture)
    res_near <- find_markers(
        exprs_sib, sib_clusters, method = "fc", top = 5,
        against = "nearest", format = "db"
    )
    res_b <- find_markers(
        exprs_sib, sib_clusters, method = "fc", top = 5,
        against = "B", format = "db"
    )
    # nearest behaves exactly like passing the resolved type explicitly
    expect_identical(
        res_near$geneSymbolmore1[res_near$cellName == "A"],
        res_b$geneSymbolmore1[res_b$cellName == "A"]
    )
    # ... and differs from using the distinct type C as the reference
    res_c <- find_markers(
        exprs_sib, sib_clusters, method = "fc", top = 5,
        against = "C", format = "db"
    )
    expect_false(identical(
        res_near$geneSymbolmore1[res_near$cellName == "A"],
        res_c$geneSymbolmore1[res_c$cellName == "A"]
    ))
})

test_that("find_markers() fc top = c(n_pos, n_neg) budgets each direction", {
    # the guard drops pan-lineage candidates from the !pos_only positive
    # pool on the main fixture, warning once about the count per type.
    # collect the warning while keeping the return value (expect_warning
    # returns the condition itself under testthat edition 3)
    w <- character(0)
    res <- withCallingHandlers(
        find_markers(
            exprs, clusters, method = "fc", top = c(5, 3),
            pos_only = FALSE, format = "db"
        ),
        warning = function(cond) {
            w <<- c(w, conditionMessage(cond))
            invokeRestart("muffleWarning")
        }
    )
    expect_match(w, "max_pct_out")
    for (ct in names(known)) {
        pos <- strsplit(
            res$geneSymbolmore1[res$cellName == ct], ","
        )[[1]]
        neg <- strsplit(
            res$geneSymbolmore2[res$cellName == ct], ","
        )[[1]]
        expect_length(pos, 5)
        expect_length(neg, 3)
    }
    # a scalar top applies to both directions
    w2 <- character(0)
    res2 <- withCallingHandlers(
        find_markers(
            exprs, clusters, method = "fc", top = 4,
            pos_only = FALSE, format = "db"
        ),
        warning = function(cond) {
            w2 <<- c(w2, conditionMessage(cond))
            invokeRestart("muffleWarning")
        }
    )
    expect_match(w2, "max_pct_out")
    expect_length(
        strsplit(res2$geneSymbolmore1[res2$cellName == "Tcell"], ",")[[1]],
        4
    )
    expect_length(
        strsplit(res2$geneSymbolmore2[res2$cellName == "Tcell"], ",")[[1]],
        4
    )
})

test_that("find_markers() fc max_pct_out guards pan-lineage genes", {
    # mA/mB/mC are type-specific markers; hA/hB/hC are over-expressed in
    # their own type but expressed in >75% of all the other cells (5 in
    # half of the cells of each other type)
    guard_counts <- matrix(
        0, 6, 80,
        dimnames = list(
            c("mA", "mB", "mC", "hA", "hB", "hC"),
            names(sib_clusters)
        )
    )
    guard_counts["mA", sib_clusters == "A"] <- 10
    guard_counts["mB", sib_clusters == "B"] <- 10
    guard_counts["mC", sib_clusters == "C"] <- 10
    guard_counts["hA", sib_clusters == "A"] <- 10
    guard_counts["hA", sib_clusters != "A"] <- 5
    guard_counts["hB", sib_clusters == "B"] <- 10
    guard_counts["hB", sib_clusters != "B"] <- 5
    guard_counts["hC", sib_clusters == "C"] <- 10
    guard_counts["hC", sib_clusters != "C"] <- 5
    exprs_g <- log1p(guard_counts)
    w <- character(0)
    res <- withCallingHandlers(
        find_markers(
            exprs_g, sib_clusters, method = "fc", top = 10, format = "db"
        ),
        warning = function(cond) {
            w <<- c(w, conditionMessage(cond))
            invokeRestart("muffleWarning")
        }
    )
    # a single aggregated warning reports one dropped candidate per type
    expect_length(w, 1)
    expect_match(w, "max_pct_out")
    expect_match(w, "A: 1, B: 1, C: 1")
    # only the type-specific markers survive
    expect_identical(res$geneSymbolmore1[res$cellName == "A"], "mA")
    expect_identical(res$geneSymbolmore1[res$cellName == "B"], "mB")
    expect_identical(res$geneSymbolmore1[res$cellName == "C"], "mC")
    # max_pct_out = 1 disables the guard: no warning, hybrids are kept
    w2 <- character(0)
    res2 <- withCallingHandlers(
        find_markers(
            exprs_g, sib_clusters, method = "fc", top = 10,
            max_pct_out = 1, format = "db"
        ),
        warning = function(cond) {
            w2 <<- c(w2, conditionMessage(cond))
            invokeRestart("muffleWarning")
        }
    )
    expect_identical(w2, character(0))
    expect_true("hA" %in% strsplit(
        res2$geneSymbolmore1[res2$cellName == "A"], ","
    )[[1]])
})

test_that("find_markers() fc `against` finds sibling-specific negatives", {
    # b1 is expressed in the sibling type B only, c1 in the distinct
    # type C only; both are low in A
    neg_counts <- matrix(
        0, 2, 80,
        dimnames = list(c("b1", "c1"), names(sib_clusters))
    )
    neg_counts["b1", sib_clusters == "B"] <- 10
    neg_counts["c1", sib_clusters == "C"] <- 10
    exprs_neg <- log1p(neg_counts)
    res_null <- find_markers(
        exprs_neg, sib_clusters, method = "fc", top = c(2, 5),
        pos_only = FALSE, format = "db"
    )
    res_ag <- find_markers(
        exprs_neg, sib_clusters, method = "fc", top = c(2, 5),
        pos_only = FALSE, against = "B", format = "db"
    )
    a_null <- strsplit(
        res_null$geneSymbolmore2[res_null$cellName == "A"], ","
    )[[1]]
    a_ag <- strsplit(
        res_ag$geneSymbolmore2[res_ag$cellName == "A"], ","
    )[[1]]
    # both genes are negative markers of A vs all the other cells ...
    expect_identical(a_null, c("b1", "c1"))
    # ... but with `against = "B"` only b1 is also low in A relative to
    # the sibling type (c1 is not expressed in B at all)
    expect_identical(a_ag, "b1")
    expect_false("c1" %in% a_ag)
})

test_that("find_markers() validates `against`, `top` and `max_pct_out`", {
    expect_error(find_markers(exprs, clusters, top = c(1, 2, 3)), "top")
    expect_error(find_markers(exprs, clusters, top = c(0, 5)), "top")
    expect_error(find_markers(exprs, clusters, top = 1.5), "top")
    expect_error(
        find_markers(exprs, clusters, max_pct_out = 0), "max_pct_out"
    )
    expect_error(
        find_markers(exprs, clusters, max_pct_out = 1.5), "max_pct_out"
    )
    expect_error(find_markers(exprs, clusters, against = 1), "against")
    expect_error(
        find_markers(exprs, clusters, against = "Rare"), "not present"
    )
    expect_error(
        find_markers(exprs, clusters, method = "presto", against = "Bcell"),
        "presto"
    )
    # a type cannot be compared against itself when it is the only type
    single <- sib_clusters[sib_clusters == "A"]
    expect_error(
        find_markers(
            exprs_sib[, names(single)], single, method = "fc",
            against = "A"
        ),
        "itself"
    )
    expect_error(
        find_markers(
            exprs_sib[, names(single)], single, method = "fc",
            against = "nearest"
        ),
        "nearest"
    )
})
