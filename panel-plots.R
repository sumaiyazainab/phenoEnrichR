# ============================================================
# panel-plots.R
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

    mapped <- enrichment_result$genes_used

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
