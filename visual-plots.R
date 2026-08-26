# ============================================================
# visualisation-plots.R
# Genomics England PanelApp visualisations
# ============================================================

# ============================================================
# Plot GEL status distribution
# ============================================================

#' Plot Genomics England PanelApp evidence status
#'
#' Creates a bar plot showing the number of Green, Amber and
#' Red genes within a PanelApp gene panel.
#'
#' @param panel_data PanelApp data frame.
#' @param status_col Column containing GEL status values.
#'
#' @return ggplot object.
#' @export

plot_panel_status <- function(panel_data,
                              status_col = "GEL_status") {

  if (!status_col %in% names(panel_data)) {
    stop("'", status_col, "' column not found.")
  }

  status <- as.character(panel_data[[status_col]])

  status[status == "3"] <- "Green"
  status[status == "2"] <- "Amber"
  status[status == "1"] <- "Red"

  df <- data.frame(
    Status = factor(
      status,
      levels = c("Green", "Amber", "Red")
    )
  )

  counts <- as.data.frame(table(df$Status))
  names(counts) <- c("Status", "Count")

  ggplot2::ggplot(
    counts,
    ggplot2::aes(
      x = Status,
      y = Count,
      fill = Status
    )
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = Count),
      vjust = -0.3,
      size = 5
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        Green = "#4CAF50",
        Amber = "#FFC107",
        Red = "#E53935"
      )
    ) +
    ggplot2::labs(
      title = "Genomics England Panel Evidence",
      x = "Evidence Category",
      y = "Number of Genes"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "none"
    )
}

# ============================================================
# Plot mapping success by GEL status
# ============================================================

#' Plot mapping success of PanelApp genes
#'
#' Shows how many genes in each GEL evidence category were
#' successfully mapped and used during enrichment analysis.
#'
#' @param panel_data Original PanelApp data frame.
#' @param enrichment_result Output from run_pheno_enrichment().
#' @param gene_col Column containing gene symbols.
#' @param status_col Column containing GEL status.
#'
#' @return ggplot object.
#' @export

plot_panel_mapping <- function(panel_data,
                               enrichment_result,
                               gene_col = "Gene Symbol",
                               status_col = "GEL_status") {

  if (!inherits(enrichment_result, "pheno_enrichment")) {
    stop("Input must be a pheno_enrichment object.")
  }

  if (!gene_col %in% names(panel_data)) {
    stop("'", gene_col, "' column not found.")
  }

  if (!status_col %in% names(panel_data)) {
    stop("'", status_col, "' column not found.")
  }

  mapped <- enrichment_result$input_gene

  panel <- panel_data

  panel$Mapping <- ifelse(
    panel[[gene_col]] %in% mapped,
    "Mapped",
    "Not mapped"
  )

  panel[[status_col]] <- factor(
    panel[[status_col]],
    levels = c("3", "2", "1"),
    labels = c("Green", "Amber", "Red")
  )

  counts <- as.data.frame(
    table(
      panel[[status_col]],
      panel$Mapping
    )
  )

  names(counts) <- c(
    "Status",
    "Mapping",
    "Count"
  )

  ggplot2::ggplot(
    counts,
    ggplot2::aes(
      x = Status,
      y = Count,
      fill = Mapping
    )
  ) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::geom_text(
      ggplot2::aes(label = Count),
      position = ggplot2::position_stack(vjust = 0.5),
      colour = "white",
      size = 4
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        Mapped = "#2E7D32",
        `Not mapped` = "grey70"
      )
    ) +
    ggplot2::labs(
      title = "Panel Gene Mapping Summary",
      subtitle = "Genes successfully used during enrichment analysis",
      x = "GEL Evidence Category",
      y = "Number of Genes",
      fill = ""
    ) +
    ggplot2::theme_minimal(base_size = 13)
}

# ============================================================
# comparison-plots.R
# Comparative visualisations for phenotype enrichment
# ============================================================


# ============================================================
# Compare enrichment methods
# ============================================================

#' Compare phenotype enrichment methods
#'
#' Compares the number of statistically significant phenotype
#' terms identified using Fisher, Parent-child and Elim
#' enrichment approaches.
#'
#' @param fisher_result pheno_enrichment object produced using Fisher.
#' @param parentchild_result pheno_enrichment object produced using Parent-child.
#' @param elim_result pheno_enrichment object produced using Elim.
#' @param fdr_cutoff Adjusted p-value significance threshold.
#' @param title Optional plot title.
#'
#' @return ggplot object.
#' @export

plot_method_comparison <- function(
    fisher_result,
    parentchild_result,
    elim_result,
    fdr_cutoff = 0.05,
    title = "Comparison of Phenotype Enrichment Methods") {

  results <- list(
    Fisher = fisher_result,
    `Parent-child` = parentchild_result,
    Elim = elim_result
  )

  valid <- vapply(
    results,
    function(x) inherits(x, "pheno_enrichment"),
    logical(1)
  )

  if (!all(valid)) {
    stop("All inputs must be pheno_enrichment objects.")
  }

  comparison <- data.frame(
    Method = names(results),
    Significant_terms = vapply(
      results,
      function(x) {
        sum(
          x$results_named$padj <= fdr_cutoff,
          na.rm = TRUE
        )
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  comparison$Method <- factor(
    comparison$Method,
    levels = c("Fisher", "Parent-child", "Elim")
  )

  ggplot2::ggplot(
    comparison,
    ggplot2::aes(
      x = Method,
      y = Significant_terms,
      fill = Method
    )
  ) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(label = Significant_terms),
      vjust = -0.3,
      size = 5
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "Number of significant phenotype terms identified ",
        "at FDR < ", fdr_cutoff
      ),
      x = "Enrichment method",
      y = "Significant phenotype terms"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold")
    )
}

# ============================================================
# Compare background gene universes
# ============================================================

#' Compare enrichment using different background universes
#'
#' Shows how the number of significant phenotype terms changes
#' when different reference gene universes are used.
#'
#' @param results Named list of pheno_enrichment objects.
#' @param fdr_cutoff Adjusted p-value significance threshold.
#' @param title Optional plot title.
#'
#' @return ggplot object.
#' @export

plot_background_comparison <- function(
    results,
    fdr_cutoff = 0.05,
    title = "Effect of Background Gene Universe on Enrichment") {

  if (is.null(names(results)) ||
      any(names(results) == "")) {
    stop("results must be a named list.")
  }

  valid <- vapply(
    results,
    function(x) inherits(x, "pheno_enrichment"),
    logical(1)
  )

  if (!all(valid)) {
    stop("All list elements must be pheno_enrichment objects.")
  }

  comparison <- data.frame(
    Background = names(results),

    Significant_terms = vapply(
      results,
      function(x) {
        sum(
          x$results_named$padj <= fdr_cutoff,
          na.rm = TRUE
        )
      },
      numeric(1)
    ),

    Universe_size = vapply(
      results,
      function(x) x$universe_size,
      numeric(1)
    ),

    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(
    comparison,
    ggplot2::aes(
      x = Background,
      y = Significant_terms,
      fill = Background
    )
  ) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(label = Significant_terms),
      vjust = -0.3,
      size = 5
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "Effect of reference gene-set selection on enrichment ",
        "(FDR < ", fdr_cutoff, ")"
      ),
      x = "Background gene universe",
      y = "Significant phenotype terms"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold")
    )
}

# ============================================================
# Compare HPO and MP enrichment
# ============================================================

#' Compare HPO and MP enrichment results
#'
#' Compares the number of significant phenotype terms and
#' input genes retained between HPO and MP enrichment.
#'
#' @param hpo_result HPO pheno_enrichment object.
#' @param mp_result MP pheno_enrichment object.
#' @param fdr_cutoff Adjusted p-value significance threshold.
#' @param title Optional plot title.
#'
#' @return ggplot object.
#' @export

plot_ontology_comparison <- function(
    hpo_result,
    mp_result,
    fdr_cutoff = 0.05,
    title = "Comparison of Human and Mouse Phenotype Enrichment") {

  if (!inherits(hpo_result, "pheno_enrichment") ||
      !inherits(mp_result, "pheno_enrichment")) {
    stop("Both inputs must be pheno_enrichment objects.")
  }

  comparison <- data.frame(
    Ontology = c("HPO", "MP"),

    Significant_terms = c(
      sum(
        hpo_result$results_named$padj <= fdr_cutoff,
        na.rm = TRUE
      ),
      sum(
        mp_result$results_named$padj <= fdr_cutoff,
        na.rm = TRUE
      )
    ),

    Genes_used = c(
      hpo_result$n_input_used,
      mp_result$n_input_used
    ),

    Universe_size = c(
      hpo_result$universe_size,
      mp_result$universe_size
    )
  )

  ggplot2::ggplot(
    comparison,
    ggplot2::aes(
      x = Ontology,
      y = Significant_terms,
      fill = Ontology
    )
  ) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = Significant_terms),
      vjust = -0.3,
      size = 5
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "Significant phenotype terms identified from HPO and MP ",
        "annotations at FDR < ", fdr_cutoff
      ),
      x = "Phenotype ontology",
      y = "Significant phenotype terms"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold")
    )
}


# ============================================================
# Top HPO vs MP enrichment terms
# ============================================================

#' Compare top HPO and MP enrichment terms
#'
#' Displays the most significantly enriched HPO and MP terms
#' side-by-side using -log10 adjusted p-values.
#'
#' @param hpo_result HPO pheno_enrichment object.
#' @param mp_result MP pheno_enrichment object.
#' @param n Number of top terms to display.
#' @param title Optional figure title.
#'
#' @return patchwork ggplot object.
#' @export

plot_hpo_mp_top_terms <- function(
    hpo_result,
    mp_result,
    n = 20,
    title = "Top HPO and MP Enriched Phenotype Terms") {

  if (!inherits(hpo_result, "pheno_enrichment") ||
      !inherits(mp_result, "pheno_enrichment")) {
    stop("Both inputs must be pheno_enrichment objects.")
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install 'patchwork' using install.packages('patchwork').")
  }

  prepare_terms <- function(result, n) {

    df <- result$results_named

    df <- df[
      !is.na(df$padj) &
        !is.na(df$term_name),
      ,
      drop = FALSE
    ]

    df$score <- -log10(
      pmax(df$padj, .Machine$double.xmin)
    )

    df <- df[
      order(df$score, decreasing = TRUE),
      ,
      drop = FALSE
    ]

    head(df, n)
  }

  hpo <- prepare_terms(hpo_result, n)
  mp <- prepare_terms(mp_result, n)

  p_hpo <- ggplot2::ggplot(
    hpo,
    ggplot2::aes(
      x = stats::reorder(term_name, score),
      y = score
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Human Phenotype Ontology (HPO)",
      x = NULL,
      y = "-log10 adjusted p-value"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  p_mp <- ggplot2::ggplot(
    mp,
    ggplot2::aes(
      x = stats::reorder(term_name, score),
      y = score
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Mammalian Phenotype Ontology (MP)",
      x = NULL,
      y = "-log10 adjusted p-value"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  p_hpo + p_mp +
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
      title = title
    )
}

# ============================================================
# Human-to-mouse orthologue mapping summary
# ============================================================

#' Summarise human-to-mouse orthologue mapping
#'
#' Calculates the number and percentage of genes from multiple
#' gene panels that map to available human-mouse orthologues.
#'
#' @param panels Named list of gene-panel data frames.
#' @param orthologs Output from read_orthologs().
#' @param gene_col Column containing human gene symbols.
#'
#' @return Mapping summary data frame.
#' @export

summarise_orthologue_mapping <- function(
    panels,
    orthologs,
    gene_col = "Gene Symbol") {

  if (is.null(names(panels))) {
    stop("panels must be a named list.")
  }

  if (!"human_gene" %in% names(orthologs)) {
    stop("orthologs must contain a human_gene column.")
  }

  available <- unique(
    toupper(trimws(orthologs$human_gene))
  )

  rows <- lapply(names(panels), function(panel_name) {

    panel <- panels[[panel_name]]

    if (!gene_col %in% names(panel)) {
      stop("'", gene_col, "' not found in ", panel_name)
    }

    genes <- unique(
      toupper(trimws(as.character(panel[[gene_col]])))
    )

    genes <- genes[
      !is.na(genes) &
        nzchar(genes)
    ]

    mapped <- sum(genes %in% available)
    total <- length(genes)
    unmapped <- total - mapped

    data.frame(
      Disease = panel_name,
      Human_genes = total,
      Successfully_mapped = mapped,
      Not_mapped = unmapped,
      Mapping_percent = round(
        100 * mapped / total,
        1
      ),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}


#' Plot human-to-mouse orthologue mapping
#'
#' @param mapping_summary Output from summarise_orthologue_mapping().
#' @param title Optional figure title.
#'
#' @return ggplot object.
#' @export

plot_orthologue_mapping_summary <- function(
    mapping_summary,
    title = "Mapping of Human Disease Genes to Mouse Orthologues") {

  mapped <- data.frame(
    Disease = mapping_summary$Disease,
    Mapping = "Mapped",
    Count = mapping_summary$Successfully_mapped
  )

  unmapped <- data.frame(
    Disease = mapping_summary$Disease,
    Mapping = "Not mapped",
    Count = mapping_summary$Not_mapped
  )

  df <- rbind(mapped, unmapped)

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Disease,
      y = Count,
      fill = Mapping
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(label = Count),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 4
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Human panel genes successfully linked to one-to-one mouse orthologues",
      x = "Disease panel",
      y = "Number of genes",
      fill = "Orthologue mapping"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold")
    )
}

# ============================================================
# HPO vs MP gene coverage
# ============================================================

#' Compare gene coverage between HPO and MP analyses
#'
#' Shows how many input genes are retained for phenotype
#' enrichment using HPO and MP annotations.
#'
#' @param hpo_results Named list of HPO pheno_enrichment objects.
#' @param mp_results Named list of MP pheno_enrichment objects.
#' @param title Optional figure title.
#'
#' @return ggplot object.
#' @export

plot_hpo_mp_gene_coverage <- function(
    hpo_results,
    mp_results,
    title = "Comparison of Gene Coverage in HPO and MP Analyses") {

  if (is.null(names(hpo_results)) ||
      is.null(names(mp_results))) {
    stop("Both result lists must be named.")
  }

  common <- intersect(
    names(hpo_results),
    names(mp_results)
  )

  if (!length(common)) {
    stop("HPO and MP result lists have no matching names.")
  }

  rows <- lapply(common, function(nm) {

    hpo <- hpo_results[[nm]]
    mp <- mp_results[[nm]]

    data.frame(
      Disease = c(nm, nm),
      Ontology = c("HPO", "MP"),
      Genes_used = c(
        hpo$n_input_used,
        mp$n_input_used
      ),
      stringsAsFactors = FALSE
    )
  })

  df <- dplyr::bind_rows(rows)

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Disease,
      y = Genes_used,
      fill = Ontology
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = Genes_used),
      position = ggplot2::position_dodge(width = 0.8),
      vjust = -0.3,
      size = 4
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Number of input genes retained in the testable HPO and MP annotation universes",
      x = "Disease panel",
      y = "Genes retained for enrichment",
      fill = "Ontology"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold")
    )
}

# ============================================================
# MP ontology DAG visualisation
# ============================================================

#' Plot a section of an ontology DAG
#'
#' Displays a selected ontology term and its ancestors to
#' demonstrate hierarchical relationships within the ontology.
#'
#' @param term_id Ontology term ID to highlight.
#' @param term_info Parsed ontology term table.
#' @param graph Child-to-parent ontology graph.
#' @param generations Number of ancestral generations to display.
#' @param title Optional plot title.
#'
#' @return ggraph object.
#' @export

plot_ontology_dag <- function(
    term_id,
    term_info,
    graph,
    generations = 4,
    title = "Example Hierarchical Structure of the MP Ontology") {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Install 'igraph' using install.packages('igraph').")
  }

  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("Install 'ggraph' using install.packages('ggraph').")
  }

  if (!term_id %in% names(graph)) {
    stop("term_id was not found in the ontology graph.")
  }

  keep <- term_id
  frontier <- term_id

  for (i in seq_len(generations)) {

    parents <- unique(
      unlist(
        graph[frontier],
        use.names = FALSE
      )
    )

    parents <- parents[
      !is.na(parents) &
        nzchar(parents)
    ]

    if (!length(parents)) {
      break
    }

    keep <- unique(c(keep, parents))
    frontier <- parents
  }

  edge_list <- lapply(keep, function(child) {

    parents <- graph[[child]]

    parents <- intersect(
      parents,
      keep
    )

    if (!length(parents)) {
      return(NULL)
    }

    data.frame(
      from = child,
      to = parents,
      stringsAsFactors = FALSE
    )
  })

  edges <- dplyr::bind_rows(edge_list)

  name_lookup <- unique(
    term_info[
      ,
      c("term_id", "term_name"),
      drop = FALSE
    ]
  )

  nodes <- data.frame(
    name = keep,
    stringsAsFactors = FALSE
  )

  nodes$label <- name_lookup$term_name[
    match(
      nodes$name,
      name_lookup$term_id
    )
  ]

  nodes$label[is.na(nodes$label)] <-
    nodes$name[is.na(nodes$label)]

  nodes$Selected <- ifelse(
    nodes$name == term_id,
    "Selected enriched term",
    "Ancestor term"
  )

  network <- igraph::graph_from_data_frame(
    edges,
    directed = TRUE,
    vertices = nodes
  )

  ggraph::ggraph(
    network,
    layout = "sugiyama"
  ) +
    ggraph::geom_edge_link(
      arrow = grid::arrow(
        length = grid::unit(3, "mm")
      ),
      end_cap = ggraph::circle(3, "mm")
    ) +
    ggraph::geom_node_label(
      ggplot2::aes(
        label = label,
        fill = Selected
      ),
      size = 3.5
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste(
        "Selected MP term:",
        term_id
      ),
      fill = ""
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      )
    )
}


# ============================================================
# Side-by-side HPO vs MP bubble plots
# ============================================================

#' Compare HPO and MP enrichment using bubble plots
#'
#' Displays HPO and MP enrichment results side-by-side.
#' The x-axis represents gene ratio, the y-axis represents
#' -log10 adjusted p-value, and bubble size represents the
#' number of input genes associated with each phenotype term.
#'
#' @param hpo_result HPO pheno_enrichment object.
#' @param mp_result MP pheno_enrichment object.
#' @param n Number of top terms to display.
#' @param label_n Number of phenotype terms to label.
#' @param title Optional overall figure title.
#'
#' @return patchwork ggplot object.
#' @export

plot_hpo_mp_bubbles <- function(
    hpo_result,
    mp_result,
    n = 20,
    label_n = 10,
    title = "Comparison of HPO and MP Phenotype Enrichment") {

  if (!inherits(hpo_result, "pheno_enrichment") ||
      !inherits(mp_result, "pheno_enrichment")) {
    stop("Both inputs must be pheno_enrichment objects.")
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package 'patchwork' is required. ",
      "Install it using install.packages('patchwork')."
    )
  }

  hpo_plot <- plot_enrichment_bubble(
    hpo_result,
    n = n,
    use_adjusted = TRUE,
    label_n = label_n
  ) +
    ggplot2::labs(
      title = "Human Phenotype Ontology (HPO)",
      subtitle = "Human phenotype enrichment",
      x = "Gene ratio",
      y = "-log10 adjusted p-value",
      size = "Input genes"
    )

  mp_plot <- plot_enrichment_bubble(
    mp_result,
    n = n,
    use_adjusted = TRUE,
    label_n = label_n
  ) +
    ggplot2::labs(
      title = "Mammalian Phenotype Ontology (MP)",
      subtitle = "Mouse phenotype enrichment via human-mouse orthologues",
      x = "Gene ratio",
      y = "-log10 adjusted p-value",
      size = "Input genes"
    )

  hpo_plot + mp_plot +
    patchwork::plot_layout(
      ncol = 2,
      guides = "collect"
    ) +
    patchwork::plot_annotation(
      title = title,
      subtitle = paste0(
        "Bubble position represents enrichment strength and gene coverage; ",
        "bubble size represents supporting input genes"
      )
    )
}

