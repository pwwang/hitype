test_that("parse_next_imm_levels() works with empty marks", {
    all_next_levels <- c("A", "B", "C")
    next_imm_levels <- parse_next_imm_levels("", all_next_levels)

    expect_equal(next_imm_levels, c(all_next_levels, UNKNOWN))
})

test_that("parse_next_imm_levels() works with NULL", {
    all_next_levels <- c("A", "B", "C")
    next_imm_levels <- parse_next_imm_levels(NULL, all_next_levels)

    expect_equal(next_imm_levels, c(all_next_levels, UNKNOWN))
})

test_that("parse_next_imm_levels() works with sole negation", {
    all_next_levels <- c("A", "B", "C")
    next_imm_levels <- parse_next_imm_levels("!", all_next_levels)

    expect_equal(next_imm_levels, "<EMPTY>")
})

test_that("parse_next_imm_levels() works with negation", {
    all_next_levels <- c("A", "B", "C")
    next_imm_levels <- parse_next_imm_levels("!A", all_next_levels)

    expect_equal(next_imm_levels, c("B", "C", UNKNOWN))
})

test_that("parse_next_imm_levels() works with positive", {
    all_next_levels <- c("A", "B", "C")
    next_imm_levels <- parse_next_imm_levels("A,B", all_next_levels)

    expect_equal(next_imm_levels, c("A", "B", UNKNOWN))
})

test_that("parse_next_levels()", {
    cell_markers <- data.frame(
        cellName = c("CD4", "CD8", "Naive", "Memory", "Activated", "Proliferating"),
        level = c(1, 1, 2, 2, 3, 3),
        nextLevels = c("", "!;!", "!Activated", "!", "", "")
    )
    max_level <- 3
    next_imm_levels <- list()
    for (i in seq_len(nrow(cell_markers))) {
        if (cell_markers$level[i] == max_level) {
            next
        }
        next_imm_levels[[cell_markers$cellName[i]]] <- parse_next_imm_levels(
            cell_markers$nextLevels[i],
            cell_markers$cellName[cell_markers$level == cell_markers$level[i] + 1]
        )
    }

    cell_name <- "Naive"
    level <- 2
    nl_marks <- "!Activated"
    next_levels <- parse_next_levels(cell_name, level, max_level, nl_marks, cell_markers, next_imm_levels)
    expect_equal(
        next_levels,
        list(Naive = c("Proliferating", UNKNOWN))
    )

    cell_name <- "CD8"
    level <- 1
    nl_marks <- "!;!"
    next_levels <- parse_next_levels(cell_name, level, max_level, nl_marks, cell_markers, next_imm_levels)
    expect_equal(
        next_levels,
        list(CD8 = "<EMPTY>")
    )

    cell_name <- "CD4"
    level <- 1
    nl_marks <- ""
    next_levels <- parse_next_levels(cell_name, level, max_level, nl_marks, cell_markers, next_imm_levels)
    expect_equal(
        next_levels,
        list(
            CD4 = list(
                Naive = c("Proliferating", UNKNOWN),
                Memory = "<EMPTY>",
                "<UNKNOWN>" = UNKNOWN
            )
        )
    )
})

test_that("gs_prepare() from data.frame", {
    gs <- gs_prepare(hitypedb_short, tissue_type = "Immune system")
    expect_equal(length(gs$gene_sets), 1)
    expect_null(gs$cell_names)
})

test_that("gs_prepare() from xlsx", {
    gs <- gs_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx")
    expect_equal(length(gs$gene_sets), 1)
    expect_null(gs$cell_names)
})

test_that("gs_prepare() from txt", {
    gs <- gs_prepare(test_path("db.txt"))
    expect_equal(length(gs$gene_sets), 3)
    expect_equal(length(gs$cell_names), 7)
})

test_that("gs_prepare() from ,; in cellName", {
    expect_error(gs_prepare(data.frame(
        cellName = "a,;b",
        geneSymbolmore1 = "a",
        geneSymbolmore2 = "a"
    )))
})

test_that("gs_prepare() no tisseType column", {
    expect_error(gs_prepare(test_path("db.txt"), tissue_type = "xyz"))
})

test_that("gs_prepare() level not from 1", {
    db <- hitypedb_short
    db$level <- 2
    expect_error(gs_prepare(db))
})

test_that("gs_prepare() level not consecutive", {
    db <- hitypedb_short
    db$level <- c(1, rep(3, nrow(db) - 1))
    expect_error(gs_prepare(db))
})