# ============================================================
# enrichment.R
# phenoEnrichR enrichment, ontology-aware methods, summaries,
# comparisons, pruning, tables, plots, and examples
# ============================================================

.make_term2genes <- function(annotations, universe, min_size = 5L) {
  ann <- validate_annotations(annotations)
  ann <- ann[ann$gene %in% universe, , drop = FALSE]
  out <- lapply(split(ann$gene, ann$term_id), unique)
  out[lengths(out) >= min_size]
}

.fisher_one <- function(selected, reference, genes_with_term) {
  selected <- unique(selected)
  reference <- setdiff(unique(reference), selected)
  a <- sum(selected %in% genes_with_term)
  b <- length(selected) - a
  c <- sum(reference %in% genes_with_term)
  d <- length(reference) - c
  mat <- matrix(c(a, b, c, d), nrow = 2L, byrow = TRUE,
                dimnames = list(Set = c("Selected", "Reference"),
                                Annotation = c("With_term", "Without_term")))
  p <- if (sum(mat) == 0L || any(rowSums(mat) == 0L)) 1 else
    stats::fisher.test(mat, alternative = "greater")$p.value
  list(a = a, b = b, c = c, d = d, p = p, odds_ratio = if ((b * c) == 0) Inf else (a * d) / (b * c))
}

.prepare_analysis_inputs <- function(genes, universe, annotations) {
  annotations <- validate_annotations(annotations)
  universe <- clean_gene_input(universe)
  supplied <- clean_gene_input(genes)
  genes_used <- intersect(supplied, universe)
  dropped <- setdiff(supplied, universe)
  if (!length(genes_used)) stop("None of the input genes were found in the selected universe.", call. = FALSE)
  if (length(dropped)) warning(length(dropped), " input gene(s) were not in the universe and were excluded.", call. = FALSE)
  list(genes = genes_used, supplied = supplied, dropped = dropped,
       universe = universe, annotations = annotations)
}

#' Classic Fisher enrichment
run_fisher_enrichment <- function(genes, universe, annotations, min_size = 5L) {
  x <- .prepare_analysis_inputs(genes, universe, annotations)
  term2genes <- .make_term2genes(x$annotations, x$universe, min_size)
  rows <- lapply(names(term2genes), function(term) {
    z <- .fisher_one(x$genes, x$universe, term2genes[[term]])
    data.frame(term_id = term,
               selected_with_term = z$a,
               selected_without_term = z$b,
               background_with_term = z$c,
               background_without_term = z$d,
               odds_ratio = z$odds_ratio,
               p_value = z$p,
               stringsAsFactors = FALSE)
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) return(out)
  out$padj <- stats::p.adjust(out$p_value, method = "BH")
  out[order(out$p_value), , drop = FALSE]
}

#' Parent-child conditional enrichment
#'
#' Tests each term within the union of genes annotated to its direct parents.
#' Root terms fall back to the full universe.
run_parentchild_enrichment <- function(genes, universe, annotations, graph, min_size = 5L) {
  x <- .prepare_analysis_inputs(genes, universe, annotations)
  term2genes <- .make_term2genes(x$annotations, x$universe, min_size)

  rows <- lapply(names(term2genes), function(term) {
    parents <- graph[[term]]
    parents <- parents[parents %in% names(term2genes)]
    conditional_universe <- if (length(parents)) {
      unique(unlist(term2genes[parents], use.names = FALSE))
    } else x$universe
    conditional_universe <- intersect(conditional_universe, x$universe)
    selected_cond <- intersect(x$genes, conditional_universe)

    if (!length(selected_cond) || length(conditional_universe) < 2L) {
      z <- list(a = 0L, b = length(selected_cond), c = 0L,
                d = length(setdiff(conditional_universe, selected_cond)),
                p = 1, odds_ratio = NA_real_)
    } else {
      z <- .fisher_one(selected_cond, conditional_universe, term2genes[[term]])
    }

    data.frame(term_id = term,
               selected_with_term = z$a,
               selected_without_term = z$b,
               background_with_term = z$c,
               background_without_term = z$d,
               conditional_universe_size = length(conditional_universe),
               n_parents_conditioned = length(parents),
               odds_ratio = z$odds_ratio,
               p_value = z$p,
               stringsAsFactors = FALSE)
  })

  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) return(out)
  out$padj <- stats::p.adjust(out$p_value, method = "BH")
  out[order(out$p_value), , drop = FALSE]
}

#' Elim-style top-down ontology enrichment
#'
#' Terms are tested from deepest to broadest. Genes belonging to a significant
#' child term are removed from its ancestors before ancestors are tested.
run_elim_enrichment <- function(genes, universe, annotations, graph,
                                min_size = 5L, elim_p_cutoff = 0.01) {
  x <- .prepare_analysis_inputs(genes, universe, annotations)
  term2genes <- .make_term2genes(x$annotations, x$universe, min_size)
  depths <- calculate_term_depth(graph)
  order_terms <- names(sort(depths[names(term2genes)], decreasing = TRUE, na.last = TRUE))
  active <- term2genes
  ancestor_map <- make_ancestor_map(order_terms, graph, include_self = FALSE)
  results <- vector("list", length(order_terms))

  for (i in seq_along(order_terms)) {
    term <- order_terms[i]
    genes_with_term <- active[[term]]
    z <- .fisher_one(x$genes, x$universe, genes_with_term)
    results[[i]] <- data.frame(
      term_id = term,
      selected_with_term = z$a,
      selected_without_term = z$b,
      background_with_term = z$c,
      background_without_term = z$d,
      term_depth = unname(depths[term]),
      odds_ratio = z$odds_ratio,
      p_value = z$p,
      stringsAsFactors = FALSE
    )

    if (is.finite(z$p) && z$p < elim_p_cutoff && length(genes_with_term)) {
      ancestors <- intersect(ancestor_map[[term]], names(active))
      for (ancestor in ancestors) active[[ancestor]] <- setdiff(active[[ancestor]], genes_with_term)
    }
  }

  out <- dplyr::bind_rows(results)
  if (!nrow(out)) return(out)
  out$padj <- stats::p.adjust(out$p_value, method = "BH")
  out[order(out$p_value), , drop = FALSE]
}

#' Optional topGO wrapper for genuine GO identifiers
#'
#' topGO is GO-specific and cannot directly analyse HP:/MP: identifiers.
#' This helper is included for benchmarking GO control datasets only.
run_topgo_control <- function(genes, universe, gene2go,
                              ontology = c("BP", "MF", "CC"),
                              algorithm = c("weight01", "elim", "weight", "parentchild", "classic")) {
  if (!requireNamespace("topGO", quietly = TRUE)) {
    stop("Install topGO with BiocManager::install('topGO').", call. = FALSE)
  }
  ontology <- match.arg(ontology)
  algorithm <- match.arg(algorithm)
  if (!all(grepl("^GO:", unique(unlist(gene2go, use.names = FALSE))))) {
    stop("run_topgo_control() only accepts GO identifiers. HPO/MP use the equivalent methods implemented in this package.", call. = FALSE)
  }
  universe <- unique(as.character(universe))
  selected <- universe %in% genes
  names(selected) <- universe
  go_data <- methods::new("topGOdata", ontology = ontology, allGenes = as.integer(selected),
                          geneSel = function(x) x == 1L,
                          annot = topGO::annFUN.gene2GO, gene2GO = gene2go)
  res <- topGO::runTest(go_data, algorithm = algorithm, statistic = "fisher")
  data.frame(term_id = names(topGO::score(res)), p_value = unname(topGO::score(res)),
             stringsAsFactors = FALSE)
}

add_term_names <- function(results, term_info) {
  if (!nrow(results)) return(results)
  names_df <- unique(term_info[, c("term_id", "term_name"), drop = FALSE])
  merge(results, names_df, by = "term_id", all.x = TRUE, sort = FALSE)
}

#' Prune redundant enriched ontology terms
#'
#' Removes an ancestor when a significant descendant has highly overlapping
#' annotated genes and is at least as significant.
prune_redundant_terms <- function(results, annotations, graph,
                                  fdr_cutoff = 0.05,
                                  overlap_threshold = 0.8,
                                  prefer = c("specific", "significant")) {
  prefer <- match.arg(prefer)
  if (!nrow(results)) return(results)
  ann <- validate_annotations(annotations)
  t2g <- lapply(split(ann$gene, ann$term_id), unique)
  sig <- results[!is.na(results$padj) & results$padj <= fdr_cutoff, , drop = FALSE]
  if (nrow(sig) < 2L) return(results)

  ancestors <- make_ancestor_map(sig$term_id, graph, include_self = FALSE)
  remove <- character()
  for (child in sig$term_id) {
    child_genes <- t2g[[child]]
    if (!length(child_genes)) next
    parent_candidates <- intersect(ancestors[[child]], sig$term_id)
    for (parent in parent_candidates) {
      parent_genes <- t2g[[parent]]
      if (!length(parent_genes)) next
      overlap <- length(intersect(child_genes, parent_genes)) / min(length(child_genes), length(parent_genes))
      child_p <- sig$padj[match(child, sig$term_id)]
      parent_p <- sig$padj[match(parent, sig$term_id)]
      if (overlap >= overlap_threshold) {
        if (prefer == "specific" || child_p <= parent_p) remove <- c(remove, parent)
      }
    }
  }
  out <- results[!results$term_id %in% unique(remove), , drop = FALSE]
  attr(out, "pruned_terms") <- unique(remove)
  out
}

filter_results_basic <- function(results,
                                 exclude_terms = c("Mode of inheritance",
                                                   "Autosomal dominant inheritance",
                                                   "Autosomal recessive inheritance",
                                                   "mammalian phenotype"),
                                 max_background = Inf,
                                 min_depth = NULL,
                                 term_depth = NULL) {
  out <- results
  if ("term_name" %in% names(out)) out <- out[!out$term_name %in% exclude_terms, , drop = FALSE]
  if ("background_with_term" %in% names(out)) out <- out[out$background_with_term <= max_background, , drop = FALSE]
  if (!is.null(min_depth) && !is.null(term_depth)) {
    out <- out[term_depth[out$term_id] >= min_depth, , drop = FALSE]
  }
  out
}

.run_enrichment_method <- function(method, genes, universe, annotations, graph,
                                   min_size, elim_p_cutoff) {
  switch(method,
         fisher = run_fisher_enrichment(genes, universe, annotations, min_size),
         parentchild = run_parentchild_enrichment(genes, universe, annotations, graph, min_size),
         elim = run_elim_enrichment(genes, universe, annotations, graph, min_size, elim_p_cutoff),
         stop("Unsupported method: ", method, call. = FALSE))
}

#' Main phenotype enrichment function
#'
#' @param genes User gene input.
#' @param ontology HPO or MP.
#' @param method fisher, parentchild, or elim.
#' @param universe Optional background. Defaults to the bundled HGNC approved protein-coding genes.
#' @param propagate Propagate annotations before testing.
#' @param prune Apply post-test redundancy pruning.
#' @param ... Resource paths and method options.
#' @return pheno_enrichment object.
#' @export
run_pheno_enrichment <- function(
    genes,
    ontology = c("HPO", "MP"),
    method = c("fisher", "parentchild", "elim"),
    gene_col = NULL,
    universe = NULL,
    propagate = TRUE,
    prune = FALSE,
    prune_overlap = 0.8,
    min_size = 5L,
    elim_p_cutoff = 0.01,
    hpo_obo_path = NULL,
    hpo_annotation_path = NULL,
    mp_obo_path = NULL,
    ortholog_path = NULL,
    mgi_genepheno_path = NULL,
    protein_coding_path = NULL) {

  ontology <- match.arg(ontology)
  method <- match.arg(method)
  genes_supplied <- clean_gene_input(genes, gene_col = gene_col)
  obo_path <- if (ontology == "HPO") hpo_obo_path else mp_obo_path

  # Package users normally use processed ontology data stored in data/*.rda.
  # A raw OBO path is only an optional developer override.
  term_info <- load_default_ontology(ontology, obo_path = obo_path)
  graph <- build_ontology_graph(term_info)
  depth <- calculate_term_depth(graph)

  ann_full <- load_default_annotations(
    ontology = ontology,
    hpo_annotation_path = hpo_annotation_path,
    ortholog_path = ortholog_path,
    mgi_genepheno_path = mgi_genepheno_path
  )
  annotations_direct <- validate_annotations(ann_full)
  annotations_used <- if (isTRUE(propagate)) propagate_annotations(annotations_direct, graph) else annotations_direct

  universe_source <- "user supplied"
  if (is.null(universe)) {
    universe <- load_default_protein_coding_genes(protein_coding_path)
    universe_source <- "HGNC approved protein-coding genes"
  }

  # Only genes that can be tested need to remain in the statistical universe.
  # This preserves the biological choice of a protein-coding reference while
  # avoiding unannotated genes creating impossible contingency-table cells.
  universe_requested <- clean_gene_input(universe)
  universe <- intersect(universe_requested, unique(annotations_used$gene))
  genes_used <- intersect(genes_supplied, universe)
  genes_unmapped <- setdiff(genes_supplied, universe)
  if (!length(genes_used)) stop("None of the input genes were found in the annotation universe.", call. = FALSE)

  results <- .run_enrichment_method(method, genes_used, universe, annotations_used,
                                    graph, min_size, elim_p_cutoff)
  results_named <- add_term_names(results, term_info)

  if (isTRUE(prune)) {
    results_named <- prune_redundant_terms(results_named, annotations_used, graph,
                                           overlap_threshold = prune_overlap)
  }

  # For MP, retain a concise orthologue mapping summary using the processed
  # package dataset. HPO does not require an orthologue step.
  mapped_input <- NULL
  if (ontology == "MP") {
    orthologs_used <- load_default_orthologs(ortholog_path)
    mapped_input <- orthologs_used[orthologs_used$human_gene %in% genes_supplied, , drop = FALSE]
  }

  object <- list(
    ontology = ontology,
    method = method,
    annotation_source = if (ontology == "HPO") "HPO" else "MGI",
    input_genes_supplied = genes_supplied,
    input_genes = genes_used,
    unmapped_input_genes = genes_unmapped,
    mapped_input = mapped_input,
    n_mapped_orthologs = if (is.null(mapped_input)) NA_integer_ else length(unique(mapped_input$human_gene)),
    n_input_supplied = length(genes_supplied),
    n_input_used = length(genes_used),
    universe = universe,
    universe_requested = universe_requested,
    universe_source = universe_source,
    universe_size = length(universe),
    universe_requested_size = length(universe_requested),
    term_info = term_info,
    graph = graph,
    term_depth = depth,
    annotations_direct = annotations_direct,
    annotations_used = annotations_used,
    propagated = isTRUE(propagate),
    pruned = isTRUE(prune),
    results = results,
    results_named = results_named,
    parameters = list(min_size = min_size, elim_p_cutoff = elim_p_cutoff,
                      prune_overlap = prune_overlap),
    generated_at = Sys.time()
  )
  class(object) <- "pheno_enrichment"
  object
}

run_hpo_enrichment <- function(genes, ...) run_pheno_enrichment(genes, ontology = "HPO", ...)
run_mp_enrichment <- function(genes, ...) run_pheno_enrichment(genes, ontology = "MP", ...)

print.pheno_enrichment <- function(x, ...) {
  sig <- if (nrow(x$results_named)) sum(x$results_named$padj <= 0.05, na.rm = TRUE) else 0L
  cat("Phenotype enrichment\n")
  cat("  Ontology:         ", x$ontology, "\n", sep = "")
  cat("  Method:           ", x$method, "\n", sep = "")
  cat("  Annotation source:", x$annotation_source, "\n")
  cat("  Input genes:      ", x$n_input_used, "/", x$n_input_supplied, " mapped\n", sep = "")
  cat("  Universe:         ", x$universe_source, "\n", sep = "")
  cat("  Testable universe:", x$universe_size, "/", x$universe_requested_size, " genes annotated\n", sep = "")
  cat("  Terms tested:     ", nrow(x$results), "\n", sep = "")
  cat("  Significant FDR<0.05: ", sig, "\n", sep = "")
  invisible(x)
}

summary.pheno_enrichment <- function(object, fdr_cutoff = 0.05, n = 10L, ...) {
  tab <- make_results_table(object, n = n, fdr_cutoff = fdr_cutoff)
  list(
    ontology = object$ontology,
    method = object$method,
    annotation_source = object$annotation_source,
    input_supplied = object$n_input_supplied,
    input_used = object$n_input_used,
    unmapped_genes = object$unmapped_input_genes,
    universe_source = object$universe_source,
    universe_requested_size = object$universe_requested_size,
    universe_size = object$universe_size,
    terms_tested = nrow(object$results),
    significant_terms = sum(object$results_named$padj <= fdr_cutoff, na.rm = TRUE),
    top_terms = tab
  )
}

#' Clean result table
make_results_table <- function(enrichment, n = 20L, fdr_cutoff = 0.05,
                               remove_broad_terms = TRUE,
                               include_genes = TRUE) {
  results <- enrichment$results_named
  if (!nrow(results)) return(data.frame())
  if (remove_broad_terms) results <- filter_results_basic(results)
  results <- results[order(results$padj, results$p_value), , drop = FALSE]
  results <- head(results, n)
  denom <- results$selected_with_term + results$selected_without_term
  ratio <- ifelse(denom > 0, results$selected_with_term / denom, NA_real_)

  overlaps <- rep(NA_character_, nrow(results))
  if (include_genes) {
    t2g <- split(enrichment$annotations_used$gene, enrichment$annotations_used$term_id)
    overlaps <- vapply(results$term_id, function(term) {
      paste(sort(intersect(enrichment$input_genes, unique(t2g[[term]]))), collapse = ";")
    }, character(1))
  }

  data.frame(
    Ontology = enrichment$ontology,
    Method = enrichment$method,
    Annotation_source = enrichment$annotation_source,
    Phenotype = results$term_name,
    Term_ID = results$term_id,
    Input_genes = results$selected_with_term,
    Gene_ratio = round(ratio, 3),
    Background_genes = results$background_with_term,
    Odds_ratio = signif(results$odds_ratio, 3),
    P_value = signif(results$p_value, 3),
    Adjusted_P_value = signif(results$padj, 3),
    Significant_FDR_0.05 = results$padj <= fdr_cutoff,
    Overlapping_genes = overlaps,
    stringsAsFactors = FALSE
  )
}

plot_enrichment <- function(enrichment, n = 10L, use_adjusted = TRUE,
                            remove_broad_terms = TRUE) {
  results <- enrichment$results_named
  if (remove_broad_terms) results <- filter_results_basic(results)
  pcol <- if (use_adjusted) "padj" else "p_value"
  df <- results[is.finite(results[[pcol]]) & results[[pcol]] > 0, , drop = FALSE]
  if (!nrow(df)) stop("No terms are available to plot.", call. = FALSE)
  df$score <- -log10(df[[pcol]])
  df <- head(df[order(df$score, decreasing = TRUE), , drop = FALSE], n)
  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(term_name, score), y = score)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(title = paste(enrichment$ontology, enrichment$method, "enrichment"),
                  x = NULL,
                  y = if (use_adjusted) "-log10 adjusted p-value" else "-log10 p-value") +
    ggplot2::theme_minimal(base_size = 12)
}

plot_enrichment_bubble <- function(enrichment, n = 15L, use_adjusted = TRUE,
                                   remove_broad_terms = TRUE, label_n = 10L) {
  results <- enrichment$results_named
  if (remove_broad_terms) results <- filter_results_basic(results)
  pcol <- if (use_adjusted) "padj" else "p_value"
  df <- results[is.finite(results[[pcol]]) & results[[pcol]] > 0, , drop = FALSE]
  if (!nrow(df)) stop("No terms are available to plot.", call. = FALSE)
  df$neg_log_p <- -log10(df[[pcol]])
  denom <- df$selected_with_term + df$selected_without_term
  df$gene_ratio <- ifelse(denom > 0, df$selected_with_term / denom, NA_real_)
  df <- head(df[order(df$neg_log_p, decreasing = TRUE), , drop = FALSE], n)
  df$label <- ifelse(seq_len(nrow(df)) <= label_n, df$term_name, "")

  p <- ggplot2::ggplot(df, ggplot2::aes(gene_ratio, neg_log_p,
                                        size = selected_with_term, label = label)) +
    ggplot2::geom_point(alpha = 0.72) +
    ggplot2::labs(title = paste(enrichment$ontology, "phenotype enrichment bubble plot"),
                  subtitle = paste("Method:", enrichment$method),
                  x = "Gene ratio", y = "-log10 adjusted p-value", size = "Input genes") +
    ggplot2::theme_minimal(base_size = 12)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p + ggrepel::geom_text_repel(max.overlaps = Inf, size = 3)
  } else {
    p + ggplot2::geom_text(check_overlap = TRUE, vjust = -0.7, size = 3)
  }
}

#' Compare background universes
compare_backgrounds <- function(genes, backgrounds, ..., n = 20L) {
  if (is.null(names(backgrounds)) || any(!nzchar(names(backgrounds)))) {
    names(backgrounds) <- paste0("background_", seq_along(backgrounds))
  }
  out <- lapply(names(backgrounds), function(nm) {
    fit <- run_pheno_enrichment(genes = genes, universe = backgrounds[[nm]], ...)
    tab <- make_results_table(fit, n = n, include_genes = FALSE)
    tab$Background_name <- nm
    tab
  })
  dplyr::bind_rows(out)
}

example_genes <- function(example = c("cardiac", "cataract")) {
  example <- match.arg(example)
  if (example == "cardiac") {
    return(c("ACTC1", "ACTN2", "BAG3", "DES", "DMD", "DSP", "LMNA",
             "MYBPC3", "MYH7", "PKP2", "TNNI3", "TNNT2", "TTN"))
  }
  c("CRYAA", "CRYAB", "CRYBB1", "CRYBB2", "CRYGC", "GJA8", "MIP", "BFSP2", "PAX6", "MAF")
}

#' Output-format reference table
#'
#' This is illustrative only; it teaches users what columns mean without
#' pretending to be a current biological result.
example_output_table <- function() {
  data.frame(
    Ontology = "HPO", Method = "parentchild", Annotation_source = "HPO",
    Phenotype = c("Example phenotype A", "Example phenotype B"),
    Term_ID = c("HP:XXXXXXX", "HP:YYYYYYY"),
    Input_genes = c(12L, 8L), Gene_ratio = c(0.60, 0.40),
    Background_genes = c(80L, 35L), Odds_ratio = c(8.2, 5.1),
    P_value = c(1e-08, 2e-05), Adjusted_P_value = c(2e-06, 0.004),
    Significant_FDR_0.05 = TRUE,
    Overlapping_genes = c("GENE1;GENE2;...", "GENE3;GENE4;..."),
    stringsAsFactors = FALSE
  )
}

#' Generate versioned precomputed example results
build_precomputed_example <- function(output_path, genes = example_genes("cardiac"), ...) {
  fit <- run_pheno_enrichment(genes = genes, ...)
  payload <- list(
    result = fit,
    table = make_results_table(fit, n = 20),
    metadata = list(generated_at = Sys.time(),
                    package_version = tryCatch(as.character(utils::packageVersion("phenoEnrichR")),
                                               error = function(e) NA_character_),
                    parameters = fit$parameters,
                    ontology = fit$ontology,
                    method = fit$method,
                    annotation_source = fit$annotation_source)
  )
  saveRDS(payload, output_path)
  invisible(payload)
}


