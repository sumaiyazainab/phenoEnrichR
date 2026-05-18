#' Run phenotype enrichment with Fisher's exact test
#'
#' Tests whether phenotype terms are overrepresented in a gene set.
#'
#' @param genes Character vector of genes of interest.
#' @param universe Character vector of background genes.
#' @param annotations Data frame with columns gene and term_id.
#' @param min_size Minimum number of genes annotated to a term.
#' @return A data frame of enrichment results.
#' @export
run_fisher_enrichment <- function(genes, universe, annotations, min_size = 1) {
  genes <- unique(genes)
  universe <- unique(universe)

  annotations <- annotations[annotations$gene %in% universe, ]

  term2genes <- split(annotations$gene, annotations$term_id)
  term2genes <- term2genes[lengths(term2genes) >= min_size]

  results <- lapply(names(term2genes), function(term) {
    genes_with_term <- unique(term2genes[[term]])

    a <- sum(genes %in% genes_with_term)
    b <- sum(!(genes %in% genes_with_term))
    c <- sum((universe %in% genes_with_term) & !(universe %in% genes))
    d <- sum(!(universe %in% genes_with_term) & !(universe %in% genes))

    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)

    p_value <- fisher.test(mat)$p.value

    data.frame(
      term_id = term,
      selected_with_term = a,
      selected_without_term = b,
      background_with_term = c,
      background_without_term = d,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })

  results <- dplyr::bind_rows(results)
  results$padj <- p.adjust(results$p_value, method = "BH")
  results <- results[order(results$p_value), ]
  rownames(results) <- NULL

  results
}



#' Basic filtering of enrichment results
#'
#' Removes overly broad or unwanted HPO terms.
#'
#' @param results Data frame of enrichment results
#' @param exclude_terms Character vector of term names to exclude
#' @param max_background Maximum allowed background_with_term count
#' @return Filtered data frame
#' @export
filter_results_basic <- function(results,
                                 exclude_terms = c(
                                   "Mode of inheritance",
                                   "Autosomal dominant inheritance",
                                   "Autosomal recessive inheritance"
                                 ),
                                 max_background = Inf) {

  out <- results

  if ("term_name" %in% names(out)) {
    out <- out[!(out$term_name %in% exclude_terms), ]
  }

  out <- out[out$background_with_term <= max_background, ]
  rownames(out) <- NULL

  out
}
