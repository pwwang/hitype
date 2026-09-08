test_that("RunHitype() only works with Seurat object", {
    expect_error(RunHitype(1))
})

test_that("RunHitype() works", {
    object <- SeuratObject::pbmc_small
    gs <- gs_prepare(
        data.frame(
            cellName = c(
                "CD4", "CD8", "Naive", "Memory", "Activated", "Proliferating"
            ),
            level = c(1, 1, 2, 2, 3, 3),
            nextLevels = c("", "!;!", "!Activated", "!", "", ""),
            geneSymbolmore1 = c(
                "CCR7,STX10,PPP1R18,ZFP36L1,GZMA",
                "TUBB1,CARD16,TMEM40,MMADHC,HLA-DMA",
                "CYB561A3,FCER1A,PRF1,HIST1H2AC,GZMM",
                "PPBP,RGS2,ARHGDIA,FGFBP2,PPP1R18",
                "HLA-DRB1,LGALS2,TAF7,CD200,MYL9",
                "NRBP1,DNAJC2,FCGRT,S100B,NCF1"
            ),
            geneSymbolmore2 = ""
        )
    )
    x <- suppressWarnings(RunHitype(object, gs))
    expect_true("CD4 Naive Proliferating" %in% x$hitype)
    expect_true("CD8" %in% x$hitype)
})
