
EMPTY <- "<EMPTY>" # nolint
UNKNOWN <- "<UNKNOWN>" # nolint

#' Split the string by the separator and trim the whitespaces.
#'
#' @param x A string.
#' @param sep The separator.
#'
#' @return A character vector.
explode <- function(x, sep = ",") {
  trimws(unlist(strsplit(x, sep, fixed = TRUE)))
}

#' Get the number of repetitive given ending characters.
#'
#' @param x A string.
#' @param ending The ending character.
#'
#' @return The number of ending character.
n_ending <- function(x, ending) {
  regex <- paste0("\\", ending, "+$")
  body <- gsub(regex, "", x)
  nchar(x) - nchar(body)
}

#' Remove the suffix of the cell name.
#'
#' @param cell_name A cell name.
#'
#' @return The cell name without suffix.
revert_cell_name <- function(cell_name) {
  sub("\\.\\.\\d+", "", cell_name)
}

#' Faster version of do.call
#'
#' @source Gmisc
#' @author Max Gordon <max at gforge.se>
#'
#' @param what The function to call
#' @param args The arguments to pass to the function
#' @param quote Whether to quote the arguments
#' @param envir The environment to look for the function in the parent frame
#'
#' @return The same results as using do.call()
do_call <- function(what, args, quote = FALSE, envir = parent.frame()) {
  # nocov start
  if (quote) {
    args <- lapply(args, enquote)
  }

  if (is.null(names(args)) ||
    is.data.frame(args)) {
    argn <- args
    args <- list()
  } else {
    # Add all the named arguments
    argn <- lapply(names(args)[names(args) != ""], as.name)
    names(argn) <- names(args)[names(args) != ""]
    # Add the unnamed arguments
    argn <- c(argn, args[names(args) == ""])
    args <- args[names(args) != ""]
  }

  if (inherits(x = what, what = "character")) {
    if (is.character(what)) {
      fn <- strsplit(what, "[:]{2,3}")[[1]]
      what <- if (length(fn) == 1) {
        get(fn[[1]], envir = envir, mode = "function")
      } else {
        get(fn[[2]], envir = asNamespace(fn[[1]]), mode = "function")
      }
    }
    call <- as.call(c(list(what), argn))
  } else if (inherits(x = what, "function")) {
    f_name <- deparse(substitute(what))
    call <- as.call(c(list(as.name(f_name)), argn))
    args[[f_name]] <- what
  } else if (inherits(x = what, what = "name")) {
    call <- as.call(c(list(what, argn)))
  }
  eval(call, envir = args, enclos = envir)
  # nocov end
}