test_that("clean_gene_input standardises genes", {
  expect_equal(clean_gene_input(c(" lmna ", "LMNA", NA, "")), "LMNA")
})

test_that("classic Fisher returns expected columns", {
  ann <- data.frame(gene = c("A", "B", "C", "D"), term_id = c("T1", "T1", "T2", "T2"))
  x <- run_fisher_enrichment(c("A", "B"), c("A", "B", "C", "D"), ann, min_size = 1)
  expect_true(all(c("term_id", "p_value", "padj", "odds_ratio") %in% names(x)))
})

test_that("propagation adds parent terms", {
  graph <- list(T1 = "ROOT", ROOT = character())
  ann <- data.frame(gene = "A", term_id = "T1")
  out <- propagate_annotations(ann, graph)
  expect_true(all(c("T1", "ROOT") %in% out$term_id))
})

test_that("parentchild method runs", {
  graph <- list(T1 = "ROOT", T2 = "ROOT", ROOT = character())
  ann <- data.frame(gene = c("A", "B", "C", "D", "A", "B", "C", "D"),
                    term_id = c("T1", "T1", "T2", "T2", rep("ROOT", 4)))
  x <- run_parentchild_enrichment(c("A", "B"), c("A", "B", "C", "D"), ann, graph, min_size = 1)
  expect_true("conditional_universe_size" %in% names(x))
})


test_that("main function has focused HPO and MP interface", {
  args <- names(formals(run_pheno_enrichment))
  expect_false("custom_annotations" %in% args)
  expect_false("annotation_source" %in% args)
  expect_true(all(c("hpo_annotation_path", "mgi_genepheno_path") %in% args))
})
