#' Read an OBO file
#'
#' Reads an ontology file in OBO format as plain text.
#'
#' @param path Path to the OBO file.
#' @return A character vector of lines from the file.
#' @export
read_obo <- function(path) {
  lines <- readLines(path, warn = FALSE)
  return(lines)
}


#' Add term names to enrichment results
#'
#' Joins term names from ontology to enrichment results.
#'
#' @param results Data frame from enrichment.
#' @param term_info Data frame from parse_obo_terms().
#' @return Results with term names added.
#' @export
add_term_names <- function(results, term_info) {

  term_names <- unique(term_info[, c("term_id", "term_name")])

  merged <- merge(
    results,
    term_names,
    by = "term_id",
    all.x = TRUE
  )

  merged
}


#' Build HPO graph from parsed ontology
#'
#' Creates a named list mapping each HPO term to its parent terms.
#'
#' @param term_info Data frame from parse_obo_terms()
#' @return A named list where each name is a term_id and each value is a character vector of parent_ids
#' @export
build_hpo_graph <- function(term_info) {
  split(term_info$parent_id, term_info$term_id)
}
