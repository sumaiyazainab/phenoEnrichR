#' Parse OBO ontology terms
#'
#' Extracts term IDs, term names, and parent relationships from an OBO file.
#'
#' @param path Path to the OBO file.
#' @return A data frame with term_id, term_name, and parent_id.
#' @export
parse_obo_terms <- function(path) {

  lines <- read_obo(path)

  term_starts <- which(lines == "[Term]")

  results <- list()

  for (i in seq_along(term_starts)) {

    start <- term_starts[i]

    end <- if (i < length(term_starts)) {
      term_starts[i + 1] - 1
    } else {
      length(lines)
    }

    block <- lines[start:end]

    id_line <- grep("^id: ", block, value = TRUE)
    name_line <- grep("^name: ", block, value = TRUE)
    parent_lines <- grep("^is_a: ", block, value = TRUE)

    if (length(id_line) == 0 || length(name_line) == 0) next

    term_id <- sub("^id: ", "", id_line[1])
    term_name <- sub("^name: ", "", name_line[1])

    if (length(parent_lines) == 0) {

      results[[length(results) + 1]] <- data.frame(
        term_id = term_id,
        term_name = term_name,
        parent_id = NA_character_,
        stringsAsFactors = FALSE
      )

    } else {

      parent_ids <- sub("^is_a: ", "", parent_lines)
      parent_ids <- sub(" !.*$", "", parent_ids)

      results[[length(results) + 1]] <- data.frame(
        term_id = term_id,
        term_name = term_name,
        parent_id = parent_ids,
        stringsAsFactors = FALSE
      )
    }
  }

  dplyr::bind_rows(results)
}
