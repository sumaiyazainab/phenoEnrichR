# ============================================================
# preprocessing.R
# phenoEnrichR preprocessing and data-source functions
# ============================================================

.assert_file <- function(path, label = "File") {
  if (is.null(path) || length(path) != 1L || !nzchar(path) || !file.exists(path)) {
    stop(label, " does not exist: ", path, call. = FALSE)
  }
  invisible(normalizePath(path, mustWork = TRUE))
}

standardise_gene_symbols <- function(x, uppercase = TRUE) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (uppercase) x <- toupper(x)
  unique(x)
}

#' Clean user gene input
#'
#' Accepts a character vector, data frame, CSV/TSV path, or one-column matrix.
#' @param genes Gene input.
#' @param gene_col Optional gene column name or number.
#' @param uppercase Convert symbols to uppercase. Recommended for human genes.
#' @return Character vector.
#' @export
clean_gene_input <- function(genes, gene_col = NULL, uppercase = TRUE) {
  if (length(genes) == 1L && is.character(genes) && file.exists(genes)) {
    ext <- tolower(tools::file_ext(genes))
    genes <- if (ext == "csv") {
      utils::read.csv(genes, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      utils::read.delim(genes, stringsAsFactors = FALSE, check.names = FALSE)
    }
  }

  if (is.matrix(genes)) genes <- as.data.frame(genes, stringsAsFactors = FALSE)

  if (is.data.frame(genes)) {
    if (is.null(gene_col)) gene_col <- names(genes)[1]
    if (is.numeric(gene_col)) gene_col <- names(genes)[gene_col]
    if (!gene_col %in% names(genes)) {
      stop("gene_col was not found in the supplied data.", call. = FALSE)
    }
    genes <- genes[[gene_col]]
  }

  standardise_gene_symbols(genes, uppercase = uppercase)
}

#' Read an OBO ontology file
#' @param path OBO path.
#' @return Character vector.
read_obo <- function(path) {
  .assert_file(path, "OBO file")
  readLines(path, warn = FALSE, encoding = "UTF-8")
}

.parse_obo_block <- function(block) {
  get_one <- function(prefix) {
    x <- grep(paste0("^", prefix), block, value = TRUE)
    if (!length(x)) NA_character_ else sub(paste0("^", prefix), "", x[1])
  }

  id <- get_one("id: ")
  name <- get_one("name: ")
  obsolete <- identical(tolower(get_one("is_obsolete: ")), "true")
  replaced_by <- get_one("replaced_by: ")

  parents <- grep("^is_a: ", block, value = TRUE)
  parents <- sub("^is_a: ", "", parents)
  parents <- sub(" !.*$", "", parents)

  alt_ids <- grep("^alt_id: ", block, value = TRUE)
  alt_ids <- sub("^alt_id: ", "", alt_ids)

  list(
    term_id = id,
    term_name = name,
    parent_ids = unique(parents),
    alt_ids = unique(alt_ids),
    is_obsolete = obsolete,
    replaced_by = replaced_by
  )
}

#' Parse ontology terms from OBO
#'
#' Reads IDs, names, is_a parents, alternate IDs, and obsolete status.
#' @param path OBO file.
#' @param include_obsolete Keep obsolete terms.
#' @return Data frame with one row per term-parent relation.
#' @export
parse_obo_terms <- function(path, include_obsolete = FALSE) {
  lines <- read_obo(path)
  starts <- which(lines == "[Term]")
  if (!length(starts)) stop("No [Term] blocks were found in the OBO file.", call. = FALSE)

  parsed <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    end <- if (i < length(starts)) starts[i + 1L] - 1L else length(lines)
    parsed[[i]] <- .parse_obo_block(lines[starts[i]:end])
  }

  parsed <- Filter(function(x) !is.na(x$term_id) && !is.na(x$term_name), parsed)
  if (!include_obsolete) parsed <- Filter(function(x) !isTRUE(x$is_obsolete), parsed)

  rows <- lapply(parsed, function(x) {
    p <- x$parent_ids
    if (!length(p)) p <- NA_character_
    data.frame(
      term_id = x$term_id,
      term_name = x$term_name,
      parent_id = p,
      is_obsolete = x$is_obsolete,
      replaced_by = x$replaced_by,
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(rows)

  alt_rows <- lapply(parsed, function(x) {
    if (!length(x$alt_ids)) return(NULL)
    data.frame(alt_id = x$alt_ids, term_id = x$term_id, stringsAsFactors = FALSE)
  })
  attr(out, "alt_id_map") <- dplyr::bind_rows(alt_rows)
  out
}

#' Build ontology graph
#' @param term_info Output from parse_obo_terms().
#' @return Named child-to-parent list.
build_ontology_graph <- function(term_info) {
  required <- c("term_id", "parent_id")
  if (!all(required %in% names(term_info))) stop("term_info is missing required columns.", call. = FALSE)
  ids <- unique(term_info$term_id)
  graph <- split(term_info$parent_id, term_info$term_id)
  graph <- lapply(graph, function(x) unique(x[!is.na(x) & nzchar(x)]))
  missing <- setdiff(ids, names(graph))
  graph[missing] <- list(character())
  graph[ids]
}

get_ancestors <- function(term, graph, include_self = TRUE) {
  seen <- character()
  visit <- function(node) {
    if (is.na(node) || !nzchar(node) || node %in% seen) return(character())
    seen <<- c(seen, node)
    parents <- graph[[node]]
    if (is.null(parents) || !length(parents)) return(node)
    unique(c(node, unlist(lapply(parents, visit), use.names = FALSE)))
  }
  ans <- visit(term)
  if (!include_self) ans <- setdiff(ans, term)
  unique(ans)
}

make_ancestor_map <- function(terms, graph, include_self = TRUE) {
  terms <- unique(as.character(terms))
  stats::setNames(lapply(terms, get_ancestors, graph = graph, include_self = include_self), terms)
}

calculate_term_depth <- function(graph) {
  memo <- new.env(parent = emptyenv())
  depth_one <- function(term, visiting = character()) {
    if (exists(term, memo, inherits = FALSE)) return(get(term, memo, inherits = FALSE))
    if (term %in% visiting) return(NA_integer_)
    parents <- graph[[term]]
    if (is.null(parents) || !length(parents)) {
      val <- 0L
    } else {
      pd <- vapply(parents, depth_one, integer(1), visiting = c(visiting, term))
      val <- if (all(is.na(pd))) NA_integer_ else 1L + max(pd, na.rm = TRUE)
    }
    assign(term, val, memo)
    val
  }
  stats::setNames(vapply(names(graph), depth_one, integer(1)), names(graph))
}

validate_annotations <- function(annotations, uppercase_genes = TRUE) {
  annotations <- as.data.frame(annotations, stringsAsFactors = FALSE)
  if (!all(c("gene", "term_id") %in% names(annotations))) {
    stop("Annotations must contain columns named gene and term_id.", call. = FALSE)
  }
  annotations <- annotations[, c("gene", "term_id"), drop = FALSE]
  annotations$gene <- trimws(as.character(annotations$gene))
  if (uppercase_genes) annotations$gene <- toupper(annotations$gene)
  annotations$term_id <- trimws(as.character(annotations$term_id))
  annotations <- annotations[
    !is.na(annotations$gene) & nzchar(annotations$gene) &
      !is.na(annotations$term_id) & nzchar(annotations$term_id), , drop = FALSE
  ]
  unique(annotations)
}

#' Propagate gene annotations to ontology ancestors
#' @param annotations gene/term_id data frame.
#' @param graph Child-to-parent graph.
#' @return Expanded annotation table.
propagate_annotations <- function(annotations, graph) {
  annotations <- validate_annotations(annotations)
  amap <- make_ancestor_map(unique(annotations$term_id), graph, include_self = TRUE)
  expanded <- lapply(seq_len(nrow(annotations)), function(i) {
    terms <- amap[[annotations$term_id[i]]]
    if (!length(terms)) terms <- annotations$term_id[i]
    data.frame(gene = annotations$gene[i], term_id = terms, stringsAsFactors = FALSE)
  })
  unique(dplyr::bind_rows(expanded))
}

#' Read HPO gene-to-phenotype annotations
#' @param path genes_to_phenotype.txt.
#' @return Standard gene/term_id table with optional term_name.
read_hpo_gene_annotations <- function(path) {
  .assert_file(path, "HPO annotation file")
  df <- utils::read.delim(path, sep = "\t", header = TRUE, quote = "",
                          stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("gene_symbol", "hpo_id")
  if (!all(required %in% names(df))) {
    stop("HPO file must contain gene_symbol and hpo_id columns.", call. = FALSE)
  }
  out <- data.frame(
    gene = toupper(df$gene_symbol),
    term_id = df$hpo_id,
    term_name = if ("hpo_name" %in% names(df)) df$hpo_name else NA_character_,
    source = "HPO",
    stringsAsFactors = FALSE
  )
  unique(out[!is.na(out$gene) & !is.na(out$term_id), , drop = FALSE])
}

#' Read one-to-one human-mouse orthologues
#'
#' Supports a clean table or the shifted export encountered in this project.
#' @param path TSV file.
#' @param human_gene_col Optional human symbol column.
#' @param mouse_gene_col Optional mouse symbol column.
#' @param mouse_mgi_col Optional MGI ID column.
#' @param human_gene_from_rownames Use row names for human symbols.
#' @return Standard orthologue table.

read_orthologs <- function(path,
                           human_gene_col = NULL,
                           mouse_gene_col = NULL,
                           mouse_mgi_col = NULL,
                           human_gene_from_rownames = NULL) {

  .assert_file(path, "Orthologue file")

  df <- utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    row.names = 1
  )

  # Auto-detect the shifted orthologue export used in this project.
  shifted <- all(
    c(
      "Human Gene Symbol",
      "Mgi Gene Acc Id",
      "Mouse Support Count Threshold"
    ) %in% names(df)
  ) &&
    mean(
      grepl("^HGNC:", df[["Human Gene Symbol"]]),
      na.rm = TRUE
    ) > 0.5 &&
    mean(
      grepl("^MGI:", df[["Mouse Support Count Threshold"]]),
      na.rm = TRUE
    ) > 0.5

  if (is.null(human_gene_from_rownames)) {
    human_gene_from_rownames <- shifted
  }

  if (shifted) {

    # In the shifted export:
    # row names = human gene symbol
    # Human Gene Symbol column = HGNC ID
    # Mgi Gene Acc Id column = mouse gene symbol
    # Mouse Support Count Threshold = MGI ID

    human_gene <- rownames(df)
    human_id <- df[["Human Gene Symbol"]]
    mouse_gene <- df[["Mgi Gene Acc Id"]]
    mouse_mgi <- df[["Mouse Support Count Threshold"]]

  } else {

    pick <- function(requested, candidates, label) {

      if (!is.null(requested)) {

        if (!requested %in% names(df)) {
          stop(
            label,
            " column not found: ",
            requested,
            call. = FALSE
          )
        }

        return(requested)
      }

      hit <- candidates[candidates %in% names(df)]

      if (!length(hit)) {
        stop(
          "Could not identify ",
          label,
          " column. Supply it explicitly.",
          call. = FALSE
        )
      }

      hit[1]
    }

    hg <- pick(
      human_gene_col,
      c(
        "human_gene",
        "Human Gene Symbol",
        "human_gene_symbol"
      ),
      "human gene"
    )

    mg <- pick(
      mouse_gene_col,
      c(
        "mouse_gene",
        "Mouse Gene Symbol",
        "mouse_gene_symbol"
      ),
      "mouse gene"
    )

    mi <- pick(
      mouse_mgi_col,
      c(
        "mouse_mgi_id",
        "Mgi Gene Acc Id",
        "MGI ID"
      ),
      "mouse MGI ID"
    )

    human_gene <- if (isTRUE(human_gene_from_rownames)) {
      rownames(df)
    } else {
      df[[hg]]
    }

    mouse_gene <- df[[mg]]
    mouse_mgi <- df[[mi]]

    human_id <- if ("human_id" %in% names(df)) {
      df$human_id
    } else {
      NA_character_
    }
  }

  out <- data.frame(
    human_gene = toupper(
      trimws(as.character(human_gene))
    ),
    human_id = as.character(human_id),
    mouse_gene = trimws(
      as.character(mouse_gene)
    ),
    mouse_mgi_id = trimws(
      as.character(mouse_mgi)
    ),
    stringsAsFactors = FALSE
  )

  out <- out[
    !is.na(out$human_gene) &
      nzchar(out$human_gene) &
      !is.na(out$mouse_gene) &
      nzchar(out$mouse_gene) &
      grepl("^MGI:", out$mouse_mgi_id),
    ,
    drop = FALSE
  ]

  unique(out)
}

#' Read MGI gene-to-MP annotations
#' @param path MGI_GenePheno report.
#' @param mgi_col MGI ID column number.
#' @param mp_col MP ID column number.
#' @return mouse_mgi_id/term_id/source data frame.
read_mgi_genepheno <- function(path, mgi_col = 7L, mp_col = 5L) {
  .assert_file(path, "MGI gene-phenotype file")
  df <- utils::read.delim(path, sep = "\t", header = FALSE, quote = "", fill = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(df) < max(mgi_col, mp_col)) stop("Unexpected MGI report format.", call. = FALSE)
  out <- data.frame(
    mouse_mgi_id = trimws(as.character(df[[mgi_col]])),
    term_id = trimws(as.character(df[[mp_col]])),
    source = "MGI",
    stringsAsFactors = FALSE
  )
  out <- out[grepl("^MGI:", out$mouse_mgi_id) & grepl("^MP:", out$term_id), , drop = FALSE]
  unique(out)
}

#' Join mouse annotations to human orthologues
#' @param orthologs Standard orthologue table.
#' @param mouse_annotations MGI mouse phenotype annotations.
#' @return Human gene to MP term table.
prepare_mp_universe <- function(orthologs, mouse_annotations) {
  needed_o <- c("human_gene", "mouse_gene", "mouse_mgi_id")
  if (!all(needed_o %in% names(orthologs))) stop("orthologs has an invalid format.", call. = FALSE)
  if (!"term_id" %in% names(mouse_annotations)) stop("mouse_annotations must contain term_id.", call. = FALSE)

  use_mgi <- "mouse_mgi_id" %in% names(mouse_annotations) &&
    any(grepl("^MGI:", mouse_annotations$mouse_mgi_id), na.rm = TRUE)

  merged <- if (use_mgi) {
    merge(orthologs, mouse_annotations, by = "mouse_mgi_id")
  } else {
    merge(orthologs, mouse_annotations, by = "mouse_gene")
  }

  out <- data.frame(
    gene = toupper(merged$human_gene),
    term_id = merged$term_id,
    mouse_gene = if ("mouse_gene" %in% names(merged)) merged$mouse_gene else NA_character_,
    mouse_mgi_id = if ("mouse_mgi_id" %in% names(merged)) merged$mouse_mgi_id else NA_character_,
    source = if ("source" %in% names(merged)) merged$source else "MP",
    stringsAsFactors = FALSE
  )
  unique(out[!is.na(out$gene) & !is.na(out$term_id), , drop = FALSE])
}

# Backward-compatible name.
prepare_mp_universe_from_mgi <- function(orthologs, mgi_annotations) {
  prepare_mp_universe(orthologs, mgi_annotations)
}

#' Retrieve a processed package dataset
#'
#' Package users normally receive these objects through the installed
#' phenoEnrichR package. During development, the function also checks the
#' global environment so the scripts can be sourced and tested directly.
#' @param name Dataset object name.
#' @return The requested dataset.
.get_package_dataset <- function(name) {
  pkg <- "phenoEnrichR"

  # Normal installed-package route: LazyData objects are available in the
  # package namespace after library(phenoEnrichR).
  if (pkg %in% loadedNamespaces()) {
    ns <- asNamespace(pkg)
    if (exists(name, envir = ns, inherits = FALSE)) {
      return(get(name, envir = ns, inherits = FALSE))
    }
  }

  # Development route: devtools::load_all() or manually loaded .rda objects.
  if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    return(get(name, envir = .GlobalEnv, inherits = FALSE))
  }

  # Explicitly ask R's data loader as a final package-data attempt.
  tmp <- new.env(parent = emptyenv())
  suppressWarnings(utils::data(list = name, package = pkg, envir = tmp))
  if (exists(name, envir = tmp, inherits = FALSE)) {
    return(get(name, envir = tmp, inherits = FALSE))
  }

  stop(
    "Processed package dataset '", name, "' is unavailable. ",
    "The package developer must run data-raw/prepare_package_data.R ",
    "before building or installing phenoEnrichR.",
    call. = FALSE
  )
}

#' Load the processed ontology resource
#'
#' Uses package-safe processed data by default. A raw OBO path can still be
#' supplied by the package developer for validation or updating resources.
#' @param ontology HPO or MP.
#' @param obo_path Optional raw OBO file override.
#' @return Parsed ontology term table.
load_default_ontology <- function(ontology, obo_path = NULL) {
  ontology <- match.arg(toupper(ontology), c("HPO", "MP"))

  if (!is.null(obo_path)) {
    return(parse_obo_terms(obo_path))
  }

  object_name <- if (ontology == "HPO") "hpo_terms" else "mp_terms"
  out <- .get_package_dataset(object_name)

  required <- c("term_id", "term_name", "parent_id")
  if (!all(required %in% names(out))) {
    stop(object_name, " has an invalid format.", call. = FALSE)
  }
  out
}

#' Load the processed annotation resource
#'
#' HPO uses the processed official HPO gene-to-phenotype table. MP uses a
#' processed human-gene-to-MP table created by joining one-to-one orthologues
#' to MGI phenotype annotations.
#'
#' Raw file paths are optional developer overrides. Ordinary package users do
#' not need to provide any paths.
#' @param ontology HPO or MP.
#' @param hpo_annotation_path Optional raw HPO annotation override.
#' @param ortholog_path Optional raw orthologue override.
#' @param mgi_genepheno_path Optional raw MGI annotation override.
#' @return Standard gene/term_id annotation table.
load_default_annotations <- function(ontology,
                                     hpo_annotation_path = NULL,
                                     ortholog_path = NULL,
                                     mgi_genepheno_path = NULL) {
  ontology <- match.arg(toupper(ontology), c("HPO", "MP"))

  if (ontology == "HPO") {
    if (!is.null(hpo_annotation_path)) {
      return(read_hpo_gene_annotations(hpo_annotation_path))
    }
    return(.get_package_dataset("hpo_annotations"))
  }

  supplied <- c(!is.null(ortholog_path), !is.null(mgi_genepheno_path))
  if (any(supplied) && !all(supplied)) {
    stop(
      "For a raw MP data override, supply both ortholog_path and ",
      "mgi_genepheno_path.",
      call. = FALSE
    )
  }

  if (all(supplied)) {
    orthologs <- read_orthologs(ortholog_path)
    mgi_annotations <- read_mgi_genepheno(mgi_genepheno_path)
    return(prepare_mp_universe(orthologs, mgi_annotations))
  }

  .get_package_dataset("mp_annotations")
}

#' Load processed one-to-one human-mouse orthologues
#'
#' This is primarily used to report MP mapping information. A raw path can be
#' supplied during data-resource development.
#' @param path Optional raw orthologue file override.
#' @return Standard orthologue table.
load_default_orthologs <- function(path = NULL) {
  if (!is.null(path)) return(read_orthologs(path))
  .get_package_dataset("human_mouse_orthologs")
}

#' Read the HGNC protein-coding gene resource
#'
#' Reads an HGNC complete-set style table and returns approved human gene
#' symbols whose locus group is protein-coding or whose locus type is
#' "gene with protein product".
#' @param path Path to the HGNC tab-delimited file.
#' @return A unique uppercase character vector of approved gene symbols.
#' @export
read_protein_coding_genes <- function(path) {
  .assert_file(path, "HGNC protein-coding gene file")
  x <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  required <- c("symbol", "status")
  if (!all(required %in% names(x))) {
    stop(
      "The HGNC file must contain at least 'symbol' and 'status' columns.",
      call. = FALSE
    )
  }

  keep <- tolower(trimws(x$status)) == "approved"
  if ("locus_group" %in% names(x)) {
    keep <- keep & tolower(trimws(x$locus_group)) == "protein-coding gene"
  } else if ("locus_type" %in% names(x)) {
    keep <- keep & tolower(trimws(x$locus_type)) == "gene with protein product"
  } else {
    stop(
      "The HGNC file must contain either 'locus_group' or 'locus_type'.",
      call. = FALSE
    )
  }

  genes <- standardise_gene_symbols(x$symbol[keep])
  if (!length(genes)) {
    stop("No approved protein-coding genes were found in the HGNC file.", call. = FALSE)
  }
  genes
}

#' Load the processed HGNC protein-coding reference universe
#'
#' @param path Optional raw HGNC file override for development.
#' @return Character vector of approved human protein-coding gene symbols.
load_default_protein_coding_genes <- function(path = NULL) {
  if (!is.null(path)) return(read_protein_coding_genes(path))
  genes <- .get_package_dataset("protein_coding_genes")
  clean_gene_input(genes)
}
