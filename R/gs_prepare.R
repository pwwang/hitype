# GNU General Public License v3.0

#' Prepare gene sets for hitype
#'
#' @author Matt Mulvahill, Panwen Wang
#'
#' @param path_to_db_file A data frame with markers or
#'   Path to the marker gene database file, it should be
#'   a tab-delimited text or excel file with the following columns:
#'   - `tissueType`: The tissue type of the cell type to be annotated.
#'       This column is required only if `cell_type` is specified.
#'   - `cellName`: The name of the cell type to be annotated.
#'   - `nextLevels`: Possible next levels of the cell type to be annotated.
#'       Introduced by `hitype`, so that we can work with hierarchical cell.
#'       Levels are separated by `;`. Cell names at each level are separated
#'       by `,`. An exclamatory mark `!` at the beginning of a level means
#'       that the cell names at this level are mutually exclusive. If the
#'       levels are less than possible next levels, then the remaining levels
#'       are all possible next levels. See the example below.
#'   - `geneSymbolmore1`: The gene symbols of the marker genes that are
#'       expected to be expressed in the cell type to be annotated.
#'       The genes can be suffixed with one or more `+`. More `+` means
#'       higher expression level. For example, `CD3E++` means the gene
#'      `CD3E` is expected to be highly expressed in the cell type.
#'   - `geneSymbolmore2`: The gene symbols of the marker genes that are
#'       expected not to be expressed in the cell type to be annotated.
#'   - `level`: The levels of the cell names. Introduced by `hitype`, so that
#'       we can work with hierarchical cell names. Different levels of
#'       `cellName`s are predicted separately. For example, If we have `CD4`
#'       as level 1 and `Naive` as level 2, then our prediction for a cell type
#'       could be `CD4 Naive`.
#'       The levels should start from 1 and be consecutive. #'
#'
#' @param tissue_type The tissue type of the cell type to be annotated.
#'   This requires the `tissueType` column in the marker gene database file.
#'   If `tissue_type` is specified, then only the cell types in the specified
#'   tissue type will be used for annotation. If `tissue_type` is not specified,
#'   then all cell types in the marker gene database file will be used for
#'   annotation.
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
        # Allow xlsx to be compatible with sctype
        ext <- tolower(tools::file_ext(path_to_db_file))
        if (ext == "xlsx" || ext == "xls") {
            cell_markers <- openxlsx::read.xlsx(path_to_db_file)
        } else {
            cell_markers <- read.table(
                path_to_db_file,
                header = TRUE,
                sep = "\t",
                stringsAsFactors = FALSE
            )
        }
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
