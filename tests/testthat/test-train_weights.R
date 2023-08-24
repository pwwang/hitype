test_that("train_weights() works", {
    path_to_gs <- data.frame(
        cellName <- c(
            "Naive CD4+ T",
            "CD14+ Mono",
            "Memory CD4+",
            "B",
            "CD8+ T",
            "FCFR3A+ Mono",
            "NK",
            "DC",
            "Platelet"
        ),
        geneSymbolmore1 <- c(
            "IL7R,CCR7",
            "CD14,LYZ",
            "IL7R,S100A4",
            "MS4A1",
            "CD8A",
            "FCGR3A,MS4A7",
            "GNLY,NKG7",
            "FCER1A,CST3",
            "PPBP"
        ),
        geneSymbolmore2 <- rep("", 9)
    )

    suppressWarnings(SeuratData::InstallData("pbmc3k"))
    pbmc_rds <- test_path("pbmc3k.rds")
    if (file.exists(pbmc_rds)) {
        pbmc <- readRDS(pbmc_rds)
    } else {
        pbmc <- pbmc3k.SeuratData::pbmc3k
        pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
        pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
        pbmc <- Seurat::NormalizeData(pbmc)
        pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
        pbmc <- Seurat::ScaleData(pbmc, features = rownames(pbmc))
        pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(object = pbmc))
        pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10)
        pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)
        pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)
        new_cluster_ids <- path_to_gs$cellName
        names(new_cluster_ids) <- levels(pbmc)
        pbmc <- Seurat::RenameIdents(pbmc, new_cluster_ids)
    }

    db <- train_weights_nn(
        path_to_gs = path_to_gs,
        exprs = pbmc
    )
    print(db)
})