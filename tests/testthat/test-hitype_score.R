set.seed(1)

test_that("valid_cell_type() works", {
    gs <- list(
        CD4 = list(Naive = "Proliferating"),
        CD8 = "<EMPTY>"
    )
    expect_equal(
        valid_cell_type(c("CD4", "Naive", "Proliferating"), gs),
        "CD4 Naive Proliferating"
    )
    expect_equal(
        valid_cell_type(c("CD4", "Naive", "Activated"), gs),
        NA_character_
    )
    expect_equal(valid_cell_type("x", gs), NA_character_)
    expect_equal(valid_cell_type(c("CD8", "x", "y"), gs), "CD8")
})

test_that("hitype_score_level() works", {
    set.seed(1)
    gs_level <- list(
        CD4 = list(
            # CD4++++, IL2RA+++, IL2RB++, IL2RG+, IL7R
            markers = c("CD4", "IL2RA", "IL2RB", "IL2RG", "IL7R", "CD68"),
            weights = c(1, .8, .6, .4, .2, -1)
        ),
        CD8 = list(
            # CD8A++++, CD8B++++, IL2RA+++, IL2RB++, IL2RG+
            markers = c("CD8A", "CD8B", "IL2RA", "IL2RB", "IL2RG", "CD68"),
            weights = c(1, 1, .8, .6, .4, -1)
        )
    )
    exprs <- data.frame(matrix(rnorm(80), nrow=8))
    rownames(exprs) <- c(
        "CD4", "CD8A", "CD8B", "IL2RA", "IL2RB", "IL2RG", "IL7R", "CD68"
    )
    score <- hitype_score_level(exprs, gs_level)
    expect_equal(dim(score), c(2, 10))
})

test_that("hitype_score() invalid expr", {
    expect_error(hitype_score(list(), NULL), "must be a matrix")
    expect_error(hitype_score(data.frame(), NULL), "is empty")
})

test_that("hitype_score() invalid gs", {
    expect_error(hitype_score(matrix(1), NULL), "must be a list")
    expect_error(hitype_score(matrix(1), list(gene_sets = list())), "is empty")
})

test_that("hitype_score() invalid scaled", {
    expect_error(
        hitype_score(matrix(1), list(1), scaled = "a"),
        "must be a logical"
    )
})

test_that("hitype_score() warns non-existing gene", {
    gs_level1 <- list(
        CD4 = list(
            markers = c("CD4", "IL2RA", "IL2RB", "IL2RG", "aaaa", "Sox2"),
            weights = c(1, .8, .6, .4, .2, -1)
        )
    )
    exprs <- data.frame(matrix(rnorm(140), nrow=14))
    rownames(exprs) <- c(
        "CD4", "IL2RA", "IL2RB", "IL2RG", "IL7R", "CD68", "CD8A",
        "CD8B", "BCL2", "CCR7", "CD27", "B3GAT1", "GZMK", "FAS"
    )
    expect_warning(hitype_score(exprs, list(gs_level1)), "Sox2 -> SOX2")
})

test_that("hitype_score() raises error for no matching genes", {
    gs_level1 <- list(
        CD4 = list(
            markers = c("a1", "a2", "a3", "a4", "a5", "a6"),
            weights = c(1, .8, .6, .4, .2, -1)
        )
    )
    exprs <- data.frame(matrix(rnorm(140), nrow=14))
    rownames(exprs) <- c(
        "CD4", "IL2RA", "IL2RB", "IL2RG", "IL7R", "CD68", "CD8A",
        "CD8B", "BCL2", "CCR7", "CD27", "B3GAT1", "GZMK", "FAS"
    )
    expect_error(
        expect_warning(hitype_score(exprs, list(gs_level1))),
        "No markers"
    )
})

test_that("hitype_score() works", {
    gs_level1 <- list(
        CD4 = list(
            markers = c("CD4", "IL2RA", "IL2RB", "IL2RG", "IL7R", "CD68"),
            weights = c(1, .8, .6, .4, .2, -1)
        ),
        CD8 = list(
            markers = c("CD8A", "CD8B", "IL2RA", "IL2RB", "IL2RG", "CD68"),
            weights = c(1, 1, .8, .6, .4, -1)
        )
    )
    gs_level2 <- list(
        Naive = list(
            markers = c("BCL2", "CCR7", "CD27", "B3GAT1"),
            weights = c(1, .8, .6, -1)
        ),
        Ttm = list(
            markers = c("CCR7", "CD27", "B3GAT1", "GZMK", "FAS"),
            weights = c(-1, 1, -1, .8, 1)
        )
    )
    exprs <- data.frame(matrix(rnorm(140), nrow=14))
    rownames(exprs) <- c(
        "CD4", "IL2RA", "IL2RB", "IL2RG", "IL7R", "CD68", "CD8A",
        "CD8B", "BCL2", "CCR7", "CD27", "B3GAT1", "GZMK", "FAS"
    )
    scores <- hitype_score(exprs, list(gs_level1, gs_level2))
    expect_equal(length(scores), 2)
    expect_equal(dim(scores[[1]]), c(2, 10))
    expect_equal(dim(scores[[2]]), c(2, 10))
})


test_that("hitype_assign_level() works", {
    cells <- paste0("C", 1:10)
    clusters <- c(3, 3, 1, 3, 1, 2, 1, 2, 2, 2)
    names(clusters) <- cells
    hitype_scores <- t(data.frame(
        CD4 = c(1.26, -0.53, -0.63, 0.91, 1.01, 0.72, -0.6, 0.54, -0.08, 1.85),
        CD8 = c(-0.85, 0.03, -1.03, -0.98, 0.0, -0.23, -0.5, 1.55,  0.09, 1.32)
    ))
    colnames(hitype_scores) <- cells
    x <- hitype_assign_level(clusters, hitype_scores, .25)
    #     Cluster  CellType Scores NCells
    # CD4        3       CD4   1.64      3
    # CD8        3 <UNKNOWN>  -1.80      3
    # CD41       1 <UNKNOWN>  -0.22      3
    # CD81       1 <UNKNOWN>  -1.53      3
    # CD42       2       CD4   3.03      4
    # CD82       2       CD8   2.73      4
    expect_equal(x$Cluster, c(3, 3, 1, 1, 2, 2))
    expect_equal(
        x$CellType,
        c("CD4", "<UNKNOWN>", "CD4", "<UNKNOWN>", "CD4", "CD8")
    )
    expect_equal(
        x$Score,
        c(0.547, 0.149, 0.332, 0.181, 0.621, 0.595),
        tolerance = 1e-3
    )
})

test_that("hitype_assign() works with single data.frame", {
    cells <- paste0("C", 1:10)
    clusters <- c(3, 3, 1, 3, 1, 2, 1, 2, 2, 2)
    names(clusters) <- cells
    hitype_scores <- t(data.frame(
        CD4 = c(1.26, -0.53, -0.63, 0.91, 1.01, 0.72, -0.6, 0.54, -0.08, 1.85),
        CD8 = c(-0.85, 0.03, -1.03, -0.98, 0.0, -0.23, -0.5, 1.55,  0.09, 1.32),
        Treg = c(1.26, -0.53, -0.63, 0.91, 1.01, 0.72, -0.6, 0.54, -0.08, 1.85)
    ))
    colnames(hitype_scores) <- cells
    x <- hitype_assign(clusters, hitype_scores)
    x <- summary(x)
    #   Level Cluster CellType Score
    #   <dbl>   <dbl> <chr>    <dbl>
    # 1     1       1 CD4      0.332
    # 2     1       2 CD4      0.621
    # 3     1       3 CD4      0.547
    expect_equal(x$Cluster, c(1, 2, 3))
    expect_equal(x$CellType, c("CD4", "CD4", "CD4"))
    expect_equal(x$Score, c(0.332, 0.621, 0.547), tolerance = 1e-3)
})

test_that("hitype_assign() stops no gs for multi-level hitype_scores", {
    expect_error(hitype_assign(1, list()), "`gs` is required")
})

test_that("hitype_assign() works with hierachical scores", {
    cells <- paste0("C", 1:10)
    clusters <- c(3, 3, 1, 3, 1, 2, 1, 2, 2, 2)
    names(clusters) <- cells

    l1_scores <- t(data.frame(
        CD4 = c(1.26, .53, .63, .91, 1.01, .72, .6, .54, .08, 1.85),
        CD8 = c(.85, .03, -1.03, .98, .0, .23, .5, 1.55, .09, 1.32),
        Treg = c(1.98, .72, .6, .55, .09, -1.03, .98, .0, .23, .5)
    ))
    colnames(l1_scores) <- cells

    l2_scores <- t(data.frame(
        Naive = c(1.01, 2.1, .72, .54, .08, 1.85, 1.26, .53, .63, .91),
        Memory = c(.11, .1, .03, -1.03, .98, .0, .23, .5, 1.55,  .09),
        Effector = c(.88, .72, .6, .55, .09, -1.03, .98, .0, .23, .5)
    ))
    colnames(l2_scores) <- cells

    l3_scores <- t(data.frame(
        Activated = c(1.23, .99, .2, .3, .4, .5, .6, .7, .8, .9),
        Resting = c(.1, .2, .3, .4, .5, .6, .7, .8, .9, -1.0)
    ))
    colnames(l3_scores) <- cells

    x <- hitype_assign(
        clusters,
        list(l1_scores, l2_scores, l3_scores),
        gs = list(
            cell_names = list(
                CD4 = list(
                    Naive = c("Activated", "Resting"),
                    Memory = c("Activated", "Resting"),
                    Effector = c("Activated")
                ),
                CD8 = list(
                    Naive = c("Activated", "Resting"),
                    Memory = c("Activated", "Resting"),
                    Effector = c("Activated")
                ),
                Treg = list(
                    Naive = c("Activated", "Resting"),
                    Memory = c("Activated", "Resting"),
                    Effector = c("Activated")
                )
            )
        ),
        threshold = 0
    )
    x <- summary(x)
    #   Cluster CellType             Score
    # *   <dbl> <chr>                <dbl>
    # 1       1 CD4 Naive Resting    0.652
    # 2       2 CD4 Naive Activated  0.679
    # 3       3 Treg Naive Activated 0.782
    expect_equal(x$Cluster, c(1, 2, 3))
    expect_equal(
        x$CellType,
        c("CD4 Naive Resting", "CD4 Naive Activated", "Treg Naive Activated")
    )
})

# test_that("hitype_assign() on real data", {
#     # Load gene sets
#     gs <- gs_prepare(hitypedb_tcell)
#     # Load expression data
#     rdsfile <- file.path("/tmp", "pbmc3kt.rds")
#     if (!file.exists(rdsfile)) {
#         download.file(
#             "https://www.dropbox.com/scl/fi/pyizrlwuklt6g9yrgf51p/pbmc3kt.rds?rlkey=fz6t9qqjjf5n8dr08vv6rhyye&dl=1",
#             rdsfile
#         )
#     }
#     pbmc3kt <- readRDS(rdsfile)

#     # Calculate cell type scores
#     scores <- suppressWarnings(  # Ignore non-exist genes
#         hitype_score(pbmc3kt@assays$RNA@scale.data, gs, scaled = TRUE)
#     )
#     cell_types <- hitype_assign(pbmc3kt$seurat_clusters, scores, gs)
#     # Cluster                      CellType
#     # 0       0                      CD4 Th17
#     # 1       1           CD8 Naïve Activated
#     # 2       2            CD4 Tscm Inhibited
#     # 3       3                     CD8 Naïve
#     # 4       4         Double Negative Naïve
#     # 5       5  CD8 Tte Terminally Exhausted
#     # 6       6                      MAIT Tem
#     # 7       7 CD4 Naïve Precursor Exhausted
#     # 8       8      MAIT Tumor Recirculating
#     expect_equal(
#         cell_types$CellType,
#         c(
#             "CD4 Th17", "CD8 Naïve Activated", "CD4 Tscm Inhibited",
#             "CD8 Naïve", "Double Negative Naïve",
#             "CD8 Tte Terminally Exhausted", "MAIT Tem",
#             "CD4 Naïve Precursor Exhausted", "MAIT Tumor Recirculating"
#         )
#     )
# })

# test_that("hitype_assign() on our data", {
#     sobjfile = ".pipen/ImmunopipeMM/SeuratClusteringOfTCells/0/output/samples.seurat.RDS"
#     sobj = readRDS(sobjfile)

#     gs = gs_prepare(hitypedb_tcell)
#     scores = hitype_score(Seurat::GetAssayData(sobj), gs)

#     cell_types = hitype_assign(sobj$seurat_clusters, scores, gs)
#     print(cell_types)
# })