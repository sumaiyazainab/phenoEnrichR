#' Plot top enriched phenotype terms
#'
#' Creates a bar plot of the top enriched phenotype terms.
#'
#' @param results Data frame of enrichment results
#' @param n Number of top terms to plot
#' @param use_adjusted Logical; use adjusted p-values if TRUE
#' @return A ggplot object
#' @export
plot_enrichment <- function(results, n = 10, use_adjusted = FALSE) {

  plot_df <- results

  score_col <- if (use_adjusted) "padj" else "p_value"

  plot_df <- plot_df[!is.na(plot_df[[score_col]]), ]
  plot_df <- plot_df[plot_df[[score_col]] > 0, ]

  plot_df$score <- -log10(plot_df[[score_col]])
  plot_df <- plot_df[order(-plot_df$score), ]
  plot_df <- head(plot_df, n)

  if (!"term_name" %in% names(plot_df)) {
    plot_df$term_name <- plot_df$term_id
  }

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = reorder(term_name, score), y = score)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Top enriched phenotype terms",
      x = "Phenotype term",
      y = if (use_adjusted) "-log10(adjusted p-value)" else "-log10(p-value)"
    )
}


#' Read gene-to-term annotations
#'
#' Reads a tab-delimited annotation table and standardizes the column names.
#'
#' @param path Path to the annotation file.
#' @param gene_col Name of the gene column.
#' @param term_col Name of the term ID column.
#' @return A data frame with columns gene and term_id.
#' @export
read_annotations <- function(path, gene_col = "gene", term_col = "term_id") {
  df <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df <- df[, c(gene_col, term_col)]
  names(df) <- c("gene", "term_id")
  unique(df)
}
