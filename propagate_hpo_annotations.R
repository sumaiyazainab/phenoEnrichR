#' Propagate HPO annotations up the ontology
#'
#' Expands gene-term annotations so that each gene is also annotated
#' to all ancestor terms of each directly annotated HPO term.
#'
#' @param annotations Data frame with columns gene and term_id
#' @param graph A named list from build_hpo_graph()
#' @return Expanded data frame with columns gene and term_id
#' @export
propagate_hpo_annotations <- function(annotations, graph) {

  expanded_list <- lapply(seq_len(nrow(annotations)), function(i) {
    gene <- annotations$gene[i]
    term <- annotations$term_id[i]

    ancestors <- get_ancestors(term, graph)

    data.frame(
      gene = gene,
      term_id = ancestors,
      stringsAsFactors = FALSE
    )
  })

  expanded <- do.call(rbind, expanded_list)
  unique(expanded)
}
