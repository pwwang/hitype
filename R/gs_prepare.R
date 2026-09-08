# GNU General Public License v3.0

#' Prepare gene sets for hitype
#'
#' @author Matt Mulvahill, Panwen Wang
#'
#' @importFrom stats na.omit
#' @importFrom utils read.table
#'
#' @param path_to_db_file A data frame with markers or a path to a marker
#'   file. Two formats are supported and auto-detected:
#'   - The hitype/ScType db format (a data.frame or a tab-delimited text or
#'     excel file with the following columns):
#'     - `tissueType`: The tissue type of the cell type to be annotated.
#'         This column is required only if `tissue_type` is specified.
#'     - `cellName`: The name of the cell type to be annotated.
#'     - `nextLevels`: Possible next levels of the cell type to be annotated.
#'         Introduced by `hitype`, so that we can work with hierarchical cell.
#'         Levels are separated by `;`. Cell names at each level are separated
#'         by `,`. An exclamatory mark `!` at the beginning of a level means
#'         that the cell names at this level are mutually exclusive. If the
#'         levels are less than possible next levels, then the remaining levels
#'         are all possible next levels. See the example below.
#'     - `geneSymbolmore1`: The gene symbols of the marker genes that are
#'         expected to be expressed in the cell type to be annotated.
#'         The genes can be suffixed with one or more `+`. More `+` means
#'         higher expression level. For example, `CD3E++` means the gene
#'        `CD3E` is expected to be highly expressed in the cell type.
#'     - `geneSymbolmore2`: The gene symbols of the marker genes that are
#'         expected not to be expressed in the cell type to be annotated.
#'     - `level`: The levels of the cell names. Introduced by `hitype`, so that
#'         we can work with hierarchical cell names. Different levels of
#'         `cellName`s are predicted separately. For example, If we have `CD4`
#'         as level 1 and `Naive` as level 2, then our prediction for a cell type
#'         could be `CD4 Naive`.
#'         The levels should start from 1 and be consecutive.
#'   - The universal marker format (see
#'     <https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/>), a long
#'     table with one row per gene per cell type:
#'     - `cell_type` (required): the cell type.
#'     - `gene` (required): the marker gene.
#'     - `direction`: `positive` or `negative` (aliases: `pos`/`neg`/`+`/`-`,
#'       case-insensitive). Optional. Defaults to `positive`.
#'     - `weight`: a numeric weight. Optional. Defaults to 1. A positive
#'       marker gets `abs(weight)`; a negative marker gets `-abs(weight)`.
#'       If no `direction` is given, the weight is used as-is (so signed
#'       weights are accepted).
#'     - `tissue`, `species`, `cancer`: optional, for filtering by
#'       `tissue_type`.
#'     - `level`: optional. The levels should start from 1 and be consecutive.
#'     Column names are case-insensitive and aliases are supported:
#'     `celltype`/`cellType`/`Type` for `cell_type`, `marker`/`gene_symbol`
#'     for `gene`, `sign` for `direction`, and `tissueType` for `tissue`.
#'     Text files with extensions `txt`/`tsv`/`csv`/`xlsx`/`xls` are read
#'     as tables; `rds`/`qs`/`qs2` files should contain a data.frame
#'     (`qs`/`qs2` require the qs2 package to be installed).
#'
#' @param tissue_type The tissue type of the cell type to be annotated.
#'   For the db format, this requires the `tissueType` column in the marker
#'   gene database file; for the universal format, the `tissue` column
#'   (or its alias `tissueType`). If `tissue_type` is specified, then only
#'   the cell types in the specified tissue type will be used for annotation.
#'   If `tissue_type` is not specified, then all cell types in the marker
#'   gene database file will be used for annotation.
#'
#' @param weight_encoding How to encoding the weights. By default, plus (+) for a
#'   positive weight; minus (-) for a negative weight and a star (*) for zero weight.
#'   No sign indicates 1. One plus indicates 2, etc. You can use a function to encode
#'   these values. For example `function(x) x*x` to make the weights squared before
#'   they are involved in the computation.
#'
#' @examples
#' # nextLevels example
#' # If we have the following cell types:
#' #   level  cellName  nextLevels
#' #   1      CD4       Naive,Memory
#' #   2      Naive
#' #   2      Memory
#' #   3      Activated
#' #   3      Proliferating
#' # Then possible final cell names are:
#' #   CD4 Naive Activated
#' #   CD4 Naive Proliferating
#' #   CD4 Memory Activated
#' #   CD4 Memory Proliferating
#' #
#' # If the `nextLevels` of CD4 is `!Naive`, then possible final cell names are:
#' #   CD4 Memory Activated
#' #   CD4 Memory Proliferating
#' #
#' # If the `nextLevels` of CD4 is `!`, then possible final cell names are:
#' #   CD4 Activated
#' #   CD4 Proliferating
#' #
#' # If the `nextLevels` of CD4 is `!Naive;!`, then possible final cell names
#' # are:
#' #   CD4 Memory
#'
#' # A gene set in the universal marker format:
#' markers <- data.frame(
#'     cell_type = c("CD4", "CD4", "CD8", "CD8", "CD8"),
#'     gene = c("CD4", "IL7R", "CD8A", "CD8B", "GZMB"),
#'     direction = c("positive", "positive", "positive", "positive", "negative"),
#'     weight = c(2, 1, 3, 1, 1)
#' )
#' gs <- gs_prepare(markers)
#' gs$gene_sets[[1]]$CD4
#' gs$gene_sets[[1]]$CD8
#'
#' @return A list with gene_sets and next_levels. The structure looks like:
#' ```r
#'   list(
#'      gene_sets = list(
#'          # level 1
#'          list( CD4 = list(markers = c(...), weights = c(...)), ... ),
#'          # level 2
#'          list( Naive = list(markers = c(...), weights = c(...)), ... )
#'      ),
#'      # All possible final cell names
#'      cell_names = list(CD4 = list(Naive = c("Activated", "Proliferating")))
#'   )
#' ```
#' @export
gs_prepare <- function(path_to_db_file, tissue_type = NULL, weight_encoding = function(x) x) {
    if (is.data.frame(path_to_db_file)) {
        cell_markers <- path_to_db_file
    } else {
        ext <- tolower(tools::file_ext(path_to_db_file))
        if (ext == "xlsx" || ext == "xls") {
            # Allow xlsx to be compatible with sctype
            cell_markers <- openxlsx::read.xlsx(path_to_db_file)
        } else if (ext == "csv") {
            cell_markers <- utils::read.csv(
                path_to_db_file,
                stringsAsFactors = FALSE
            )
        } else if (ext == "rds") {
            cell_markers <- readRDS(path_to_db_file)
            if (!is.data.frame(cell_markers)) {
                stop(
                    "The `rds` file should contain a data.frame of markers."
                )
            }
        } else if (ext == "qs2" || ext == "qs") {
            if (!requireNamespace("qs2", quietly = TRUE)) {
                stop("The package `qs2` is required to read `.qs2` or `.qs` files.")
            }
            cell_markers <- qs2::qs_read(path_to_db_file)
            if (!is.data.frame(cell_markers)) {
                stop(
                    "The `qs2`/`qs` file should contain a data.frame of markers."
                )
            }
        } else {
            cell_markers <- read.table(
                path_to_db_file,
                header = TRUE,
                sep = "\t",
                stringsAsFactors = FALSE
            )
        }
    }

    # Detect the biopipen universal marker format (a long table with
    # `cell_type` and `gene` columns), which the hitype/ScType db columns
    # never have
    if (is_universal_df(cell_markers)) {
        return(gs_prepare_universal(
            canonicalize_marker_cols(cell_markers),
            tissue_type,
            weight_encoding
        ))
    }

    # Filter by tissue type
    if (!is.null(tissue_type) && is.null(cell_markers$tissueType)) {
        stop(
            "The marker gene db file does not have the `tissueType` column."
        )
    }
    if (!is.null(tissue_type)) {
        cell_markers <- cell_markers[cell_markers$tissueType == tissue_type, ]
    }

    # Set default level
    if (!"level" %in% names(cell_markers)) {
        cell_markers$level <- 1
    }

    # Check if level starts from 1 and is consecutive
    if (min(cell_markers$level) != 1) {
        stop("Level should start from 1.")
    }
    if (length(unique(cell_markers$level)) != max(cell_markers$level)) {
        stop("Level should be consecutive.")
    }

    if (any(grepl("[,;]", cell_markers$cellName))) {
        stop("The cell names should not contain `,` or `;`. ")
    }
    cell_markers$cellName <- paste(
        cell_markers$cellName,
        cell_markers$level,
        sep = ".."
    )
    cell_markers$geneSymbolmore1 <- as.character(cell_markers$geneSymbolmore1)
    cell_markers$geneSymbolmore2 <- as.character(cell_markers$geneSymbolmore2)

    # Get the max number of +/- in all geneSymbolmore1
    max_plus <- max(n_ending(explode(cell_markers$geneSymbolmore1), "+"))
    max_minus <- max(n_ending(explode(cell_markers$geneSymbolmore1), "-"))

    gene_sets <- lapply(
        split(cell_markers, cell_markers$level),
        function(x) {
            lapply(
                split(x, revert_cell_name(x$cellName)),
                function(y) {
                    markers1 <- na.omit(unique(explode(y$geneSymbolmore1)))
                    markers2 <- na.omit(unique(explode(y$geneSymbolmore2)))
                    markers2 <- setdiff(markers2, markers1)
                    genes <- gsub("\\++$|\\-+|\\*$", "", markers1)
                    weights <- unlist(sapply(markers1, function(x) {
                        if (endsWith(x, "-")) {
                            -n_ending(x, "-")
                        } else if (endsWith(x, "+")) {
                            # +1 to avoid 0
                            n_ending(x, "+") + 1
                        } else if (endsWith(x, "*")) {
                            0
                        } else {
                            1
                        }
                    }))
                    genes <- c(genes, markers2)
                    weights <- c(weights, rep(-1, length(markers2)))
                    weights <- weight_encoding(weights)
                    list(markers = genes, weights = unname(weights))
                }
            )
        }
    )

    if (is.null(cell_markers$nextLevels) || all(cell_markers$level == 1)) {
        return(list(gene_sets = gene_sets, cell_names = NULL))
    }
    max_level <- max(cell_markers$level)
    # Parse next immediate levels of cell types
    next_imm_levels <- list()
    for (i in seq_len(nrow(cell_markers))) {
        level <- cell_markers$level[i]
        if (level == max_level) {
            next
        }
        next_levels <- cell_markers$nextLevels[i]
        if (!(
            is.na(next_levels) ||
            is.null(next_levels) ||
            (is.character(next_levels) && nchar(next_levels) == 0)
        )) {
            next_levels <- explode(next_levels, ";")
            next_levels <- sapply(
                next_levels,
                function(x) {
                    if (x == "!") {
                        x
                    } else {
                        paste(explode(x), level + 1, sep = "..", collapse = ",")
                    }
                }
            )
            next_levels <- paste(next_levels, collapse = ";")
        }
        cell_markers$nextLevels[i] <- next_levels

        next_imm_levels[[cell_markers$cellName[i]]] <- parse_next_imm_levels(
            cell_markers$nextLevels[i],
            cell_markers$cellName[
                cell_markers$level == cell_markers$level[i] + 1
            ]
        )
    }

    level_1_cell_markers <- cell_markers[
        cell_markers$level == 1, , drop = FALSE
    ]
    cell_names <- list()
    for (i in seq_len(nrow(level_1_cell_markers))) {
        cell_name <- level_1_cell_markers$cellName[i]
        nl_marks <- explode(level_1_cell_markers$nextLevels[i], ";")
        cnames <- parse_next_levels(
            cell_name,
            1,
            max_level = max_level,
            nl_marks = nl_marks,
            cell_markers = cell_markers,
            next_imm_levels = next_imm_levels
        )
        cell_names <- c(cell_names, cnames)
    }
    list(gene_sets = gene_sets, cell_names = cell_names)
}

#' Canonical column names of the universal marker format
#'
#' Canonical names win over aliases. The hitype/ScType db columns
#' (cellName, geneSymbolmore1/2, tissueType...) never collide with these.
#' @keywords internal
UNIVERSAL_COL_ALIASES <- list( # nolint
    cell_type = c("cell_type", "celltype", "type"),
    gene = c("gene", "marker", "gene_symbol"),
    direction = c("direction", "sign"),
    weight = c("weight"),
    species = c("species"),
    cancer = c("cancer"),
    tissue = c("tissue", "tissueType"),
    level = c("level")
)

#' Rename case-insensitive alias columns to their canonical names
#' @keywords internal
#' @param df A data.frame of markers.
#' @return The data.frame with canonical column names.
canonicalize_marker_cols <- function(df) {
    cn <- colnames(df)
    low_cn <- tolower(cn)
    for (canonical in names(UNIVERSAL_COL_ALIASES)) {
        if (canonical %in% cn) {
            next
        }
        hit <- which(low_cn %in% tolower(UNIVERSAL_COL_ALIASES[[canonical]]))[1]
        if (!is.na(hit)) {
            cn[hit] <- canonical
        }
    }
    colnames(df) <- cn
    df
}

#' A frame is in the universal marker format iff it has (case-insensitively)
#' a cell_type-family and a gene-family column. db frames never do.
#' @keywords internal
#' @param df A data.frame of markers.
#' @return `TRUE` if the frame is in the universal marker format.
is_universal_df <- function(df) {
    if (!is.data.frame(df)) {
        return(FALSE)
    }
    cn <- tolower(colnames(df))
    any(cn %in% c("cell_type", "celltype", "type")) &&
        any(cn %in% c("gene", "marker", "gene_symbol"))
}

#' Prepare a gene set in the universal marker format
#'
#' The universal marker format is a long table (one row per gene per cell
#' type) with `cell_type` and `gene` columns, as used by biopipen
#' (<https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/>). The columns
#' are canonicalized by [canonicalize_marker_cols()] before this is called.
#'
#' @keywords internal
#'
#' @param cm The canonicalized marker table
#' @param tissue_type The tissue type to filter the markers by
#' @param weight_encoding The weight encoding function
#'
#' @return A list with `gene_sets` and `cell_names` (NULL)
gs_prepare_universal <- function(cm, tissue_type = NULL, weight_encoding = function(x) x) {
    # Filter by tissue type
    if (!is.null(tissue_type)) {
        if (is.null(cm$tissue)) {
            stop("The marker table does not have the `tissue` column.")
        }
        cm <- cm[cm$tissue == tissue_type, , drop = FALSE]
    }

    # Drop rows with missing/empty cell types or genes
    keep <- !is.na(cm$cell_type) & !is.na(cm$gene) &
        trimws(as.character(cm$cell_type)) != "" &
        trimws(as.character(cm$gene)) != ""
    cm <- cm[keep, , drop = FALSE]
    cm$cell_type <- trimws(as.character(cm$cell_type))
    cm$gene <- trimws(as.character(cm$gene))
    if (any(grepl("[,;]", cm$cell_type))) {
        stop("The cell names should not contain `,` or `;`. ")
    }

    # Validate levels, the same as the db format
    if (is.null(cm$level)) {
        cm$level <- 1
    }
    if (min(cm$level) != 1) {
        stop("Level should start from 1.")
    }
    if (length(unique(cm$level)) != max(cm$level)) {
        stop("Level should be consecutive.")
    }

    # Normalize directions (positive/negative, aliases: pos/neg/+/-
    # case-insensitive). NA/"" means no direction for that row.
    has_dir <- !is.null(cm$direction)
    if (has_dir) {
        direction <- tolower(trimws(as.character(cm$direction)))
        direction[is.na(direction) | direction == ""] <- NA_character_
        bad <- !is.na(direction) &
            !direction %in% c("positive", "pos", "neg", "negative", "+", "-")
        if (any(bad)) {
            stop(
                "Invalid `direction` value(s): ",
                paste(unique(direction[bad]), collapse = ", "),
                ". Accepted values: positive/negative (aliases: pos/neg/+/-)."
            )
        }
    }

    # Decode weights per row: the direction is authoritative for the sign.
    # - positive marker: abs(weight)
    # - negative marker: -abs(weight)
    # - no direction but a weight column: weight as-is (already signed)
    # - neither: 1
    has_wt <- !is.null(cm$weight)
    weight <- if (has_wt) {
        raw_wt <- as.character(cm$weight)
        w <- suppressWarnings(as.numeric(raw_wt))
        # Rows with an empty/missing weight fall back to the default of 1,
        # as if the `weight` column were absent for that row
        w[is.na(raw_wt) | trimws(raw_wt) == ""] <- 1
        if (anyNA(w)) {
            stop("The `weight` column contains non-numeric values.")
        }
        w
    } else {
        rep(1, nrow(cm))
    }
    # Rows with a valid direction: the direction is authoritative for the
    # sign (positive -> abs(weight), negative -> -abs(weight)).
    # Rows without a direction (column absent, or NA/""): signed weight
    # as-is if a weight column exists, otherwise 1.
    dir_neg <- rep(FALSE, nrow(cm))
    dir_pos <- rep(FALSE, nrow(cm))
    if (has_dir) {
        dir_neg <- !is.na(direction) &
            direction %in% c("negative", "neg", "-")
        dir_pos <- !is.na(direction) &
            direction %in% c("positive", "pos", "+")
    }
    signed <- rep(1, nrow(cm))
    if (has_wt) {
        signed <- weight
    }
    signed[dir_neg] <- -abs(weight[dir_neg])
    signed[dir_pos] <- abs(weight[dir_pos])

    rows <- data.frame(
        level = cm$level,
        cell_type = cm$cell_type,
        gene = cm$gene,
        signed = signed,
        stringsAsFactors = FALSE
    )
    gene_sets <- lapply(
        split(rows, rows$level),
        function(x) {
            lapply(
                split(x, x$cell_type, drop = TRUE),
                function(y) {
                    y <- y[!duplicated(y$gene), , drop = FALSE]
                    list(
                        markers = y$gene,
                        weights = unname(weight_encoding(y$signed))
                    )
                }
            )
        }
    )
    list(gene_sets = gene_sets, cell_names = NULL)
}

#' Parse all next levels of current cell name
#'
#' @keywords internal
#'
#' @param cell_name The current cell name
#' @param level The current level
#' @param max_level The max level
#' @param nl_marks The next level marks
#' @param cell_markers The cell markers data frame
#' @param next_imm_levels The next immediate levels
#'
#' @return A list of all possible cell names
parse_next_levels <- function(
    cell_name,
    level,
    max_level,
    nl_marks,
    cell_markers,
    next_imm_levels
) {
    nl_names <- parse_next_imm_levels(nl_marks[1], next_imm_levels[[cell_name]])
    if (EMPTY %in% nl_names) {
        return(as.list(structure(
            list(EMPTY),
            names = revert_cell_name(cell_name)
        )))
    }
    if (level == max_level - 1) {
        return(as.list(
            structure(
                list(revert_cell_name(nl_names)),
                names = revert_cell_name(cell_name)
            )
        ))
    }
    cell_names <- list()
    for (nl_name in nl_names) {
        cell_names <- c(
            cell_names,
            parse_next_levels(
                nl_name,
                level + 1,
                max_level,
                nl_marks[-1],
                cell_markers,
                next_imm_levels
            )
        )
    }
    return (as.list(structure(
        list(cell_names),
        names = revert_cell_name(cell_name)
    )))
}

#' Parse next immediate levels of cell types
#'
#' @keywords internal
#'
#' @param nl_marks The marks for the next levels of cell types
#' @param all_next_levels All possible next levels of cell types
#'
#' @return A vector of possible next immediate levels of cell types
parse_next_imm_levels <- function(nl_marks, all_next_levels) {
    if (length(all_next_levels) == 1 && all_next_levels == EMPTY) {
        return(EMPTY)
    }
    if (
        is.null(nl_marks) ||
        is.na(nl_marks) ||
        length(nl_marks) == 0 ||
        nchar(nl_marks) == 0
    ) {
        return(union(all_next_levels, UNKNOWN))
    }
    next_immediate <- explode(nl_marks, ";")[1]
    if (next_immediate == "!") {
        return(EMPTY)
    }
    if (startsWith(nl_marks, "!")) {
        negated <- TRUE
        nl_marks <- gsub("^!", "", nl_marks)
    } else {
        negated <- FALSE
    }
    nl_marks <- explode(nl_marks)
    if (!negated) {
        return(union(nl_marks, UNKNOWN))
    }
    return(union(setdiff(all_next_levels, nl_marks), UNKNOWN))
}
