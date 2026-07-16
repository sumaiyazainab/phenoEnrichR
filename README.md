# phenoEnrichR
This project will develop an R package for enrichment analysis using Human and Mouse Phenotype Ontologies (HPO/MP). Given gene sets from RNA-seq or GWAS, it will identify overrepresented phenotypes, assess ontology/background effects, generate visualisations, and potentially include a Shiny interface for reproducible analysis.

The aim of this project is to develop an R package that enables researchers to perform enrichment analysis using the Human Phenotype Ontology (HPO) and Mouse Phenotype (MP) ontologies. Given a set of genes identified through experimental or computational approaches such as RNA-seq experiments or GWAS studies, enrichment analysis of abnormal phenotypes described in mouse and humans will allow further characterisation of a gene set. This approach is similar to other widely implemented methods such as Gene Ontology and pathway enrichment analysis.
Given a set of genes as input, the package will identify overrepresented phenotypic terms, providing summary statistics and visualizations to facilitate interpretation. The tool will integrate existing ontologies and gene-phenotype associations in a user-friendly and reproducible workflow, allowing users to explore potential gene-phenotype relationships efficiently.
The project involves:
Developing the statistical framework to perform the enrichment analysis, accounting for the hierarchical structure of the ontologies.
Assessing the impact of different background sets, annotation sets, and ontology pruning approaches on the enrichment analysis.
Developing a set of visualisations to represent the results of the enrichment analysis.
Potentially deploying a Shiny-based graphical interface for the R package.

phenoEnrichR was developed to facilitate phenotype enrichment analysis using both human and mouse phenotype ontologies. The package enables researchers to investigate whether a set of genes is significantly associated with particular phenotypic abnormalities, providing biological context for experimental or computational gene lists.

The workflow supports:

HPO Enrichment Analysis
Uses Human Phenotype Ontology (HPO) annotations.
Identifies enriched human disease phenotypes associated with the input genes.
MP Enrichment Analysis
Maps human genes to high-confidence mouse orthologues.
Uses Mouse Genome Informatics (MGI) phenotype annotations.
Identifies enriched mouse phenotypes that may provide functional insight into human disease genes.
Ontology-Aware Analysis
Parses ontology structures directly from OBO files.
Propagates annotations through parent-child relationships.
Supports comparison of ontology pruning and annotation propagation strategies.
Visualisation and Reporting
Enrichment summary tables.
Bubble plots.
Bar plots.
Comparison of alternative background gene universes.


This version intentionally supports two built-in workflows only:

1. **HPO enrichment** using the HPO ontology and official HPO gene-to-phenotype annotations.
2. **MP enrichment** using the MP ontology, one-to-one human-mouse orthologues, and MGI phenotype annotations.

Users choose the ontology through `ontology = "HPO"` or `ontology = "MP"`. The package does not require IMPC files and does not accept custom annotation tables. This keeps the public workflow focused and reproducible.

## HPO example

```r
fit <- run_pheno_enrichment(
  genes = example_genes("cardiac"),
  ontology = "HPO",
  method = "parentchild",
  hpo_obo_path = "data-raw/hp.obo",
  hpo_annotation_path = "data-raw/genes_to_phenotype.txt"
)
```

## MP example

```r
fit <- run_pheno_enrichment(
  genes = example_genes("cardiac"),
  ontology = "MP",
  method = "parentchild",
  mp_obo_path = "data-raw/MPheno_OBO.ontology.txt",
  ortholog_path = "data-raw/One To one orthologs2026-06-08.tsv",
  mgi_genepheno_path = "data-raw/MGI_GenePheno.rpt.txt"
)




## Implemented phenoEnrichR


1. What preprocessing.R does
This script handles everything that must happen before enrichment testing.
It cleans the user's gene list
The function:

clean_gene_input()

allows users to provide genes as:
* a character vector;
* a data frame;
* a matrix;
* a CSV file;
* a TSV file.
For example:

genes <- c("LMNA", "MYH7", "TNNT2")

or:

genes <- read.csv("my_genes.csv")

The function then:
* removes missing values;
* removes empty entries;
* removes duplicate genes;
* trims spaces;
* converts human gene symbols to uppercase.
For example:

c(" lmna ", "LMNA", "", NA, "myh7")

becomes:

c("LMNA", "MYH7")


It reads and checks raw files
The internal function:

.assert_file()

checks that a requested file exists before trying to read it.
Instead of producing a confusing R error, it gives a clearer message such as:

OBO file does not exist: data-raw/hp.obo

This is mainly used during data preparation or when a developer supplies a newer raw file.

It parses HPO and MP ontology files
The function:

parse_obo_terms()

reads an OBO ontology file and extracts:

term_id
term_name
parent_id
is_obsolete
replaced_by

For example:

term_id       HP:0001644
term_name     Dilated cardiomyopathy
parent_id     HP:0001638

The parser also recognises:
* alternative term IDs;
* obsolete terms;
* replacement term IDs;
* multiple parents.
This is important because HPO and MP are directed acyclic graphs rather than simple flat lists.
By default, obsolete ontology terms are removed.

It builds the ontology graph
The function:

build_ontology_graph()

converts the parsed ontology into a child-to-parent structure.
For example:

HP:0001644
    ↓ parent
HP:0001638
    ↓ parent
HP:0001626

This graph is required for:
* finding ancestors;
* propagating annotations;
* parent–child testing;
* elim-style testing;
* calculating ontology depth;
* pruning redundant results.

It identifies ontology ancestors
The function:

get_ancestors()

finds all broader terms above a phenotype.
For example, a specific cardiomyopathy term may return:

Dilated cardiomyopathy
Cardiomyopathy
Abnormal heart morphology
Abnormality of the cardiovascular system
Phenotypic abnormality

The function includes protection against repeatedly visiting the same term.

It creates an ancestor lookup table
The function:

make_ancestor_map()

calculates the ancestors once and stores them.
This improves performance because the package does not need to repeatedly traverse the ontology for every gene–term annotation.

It calculates ontology depth
The function:

calculate_term_depth()

calculates how far each term is from the root of the ontology.
For example:

Phenotypic abnormality                     depth 0
Abnormal cardiovascular phenotype         depth 1
Abnormal heart morphology                 depth 2
Cardiomyopathy                             depth 3
Dilated cardiomyopathy                    depth 4

Ontology depth can help identify broad and specific terms and can support pruning or result filtering.

It validates annotation tables
The function:

validate_annotations()

ensures that annotation tables contain the two required columns:

gene
term_id

It also:
* converts values to character format;
* removes missing values;
* removes empty values;
* standardises gene symbols;
* removes duplicates.
This means the enrichment functions receive a consistent annotation format regardless of whether the source is HPO or MGI.

It propagates annotations through the ontology
The function:

propagate_annotations()

takes a direct gene–phenotype annotation and adds its ancestors.
For example:

LMNA → Dilated cardiomyopathy

may become:

LMNA → Dilated cardiomyopathy
LMNA → Cardiomyopathy
LMNA → Abnormal heart morphology
LMNA → Abnormal cardiovascular phenotype
LMNA → Phenotypic abnormality

This accounts for the hierarchical structure of HPO and MP.
It also removes duplicate gene–term pairs.

It reads HPO annotations
The function:

read_hpo_gene_annotations()

reads the HPO genes_to_phenotype.txt file.
It converts the source columns into the standard package format:

gene
term_id
term_name
source

For example:

LMNA    HP:0001644    Dilated cardiomyopathy    HPO


It reads the human–mouse orthologue table
The function:

read_orthologs()

reads the one-to-one human–mouse orthologue mapping.
It produces a consistent table containing information such as:

human_gene
mouse_gene
mouse_mgi_id

The function also contains logic for the unusual shifted orthologue file format encountered in your project.
This allows the MP workflow to begin with human genes while still using mouse phenotype annotations.

It reads MGI gene–phenotype annotations
The function:

read_mgi_genepheno()

reads the MGI gene–phenotype report and extracts:

mouse_mgi_id
term_id

It validates that:
* mouse gene IDs begin with MGI:;
* phenotype terms begin with MP:.
This helps detect unexpected file formats.

It creates human gene-to-MP annotations
The function:

prepare_mp_universe_from_mgi()

joins:

human–mouse orthologues

with:

mouse MGI–MP annotations

The result is a human gene-to-MP term table:

human gene
    ↓
mouse orthologue
    ↓
MP phenotype

For example:

LMNA → mouse Lmna → MP phenotype terms

Importantly, the human gene remains the unit used in the enrichment analysis. This avoids treating the mouse gene as a separate input gene.

It loads package-safe datasets
The most important recent improvement is the addition of functions such as:

load_default_ontology()
load_default_annotations()
load_default_orthologs()

These functions load processed datasets from the installed package.
The package expects datasets such as:

hpo_terms
hpo_annotations
mp_terms
mp_annotations
human_mouse_orthologs

Therefore, ordinary users do not need to provide:

"data-raw/hp.obo"

or:

"data-raw/MGI_GenePheno.rpt.txt"

They will simply run:

run_pheno_enrichment(
  genes = genes,
  ontology = "HPO"
)

This change makes the analysis package-safe.

2. What enrichment.R does
This script performs the statistical analysis after the ontology and annotation data have been prepared.
It prepares term-to-gene sets
The function:

.make_term2genes()

converts the annotation table from:

gene    term_id

into:

term_id → genes

For example:

HP:0001644 → LMNA, MYH7, TNNT2, BAG3

It also:
* restricts annotations to the selected universe;
* removes duplicates;
* excludes terms containing fewer genes than min_size.

It constructs Fisher contingency tables
The function:

.fisher_one()

constructs the four values required for Fisher's Exact Test:
	Has phenotype term	Does not have phenotype term
Selected genes	a	b
Other universe genes	c	d
It then calculates:
* the one-sided Fisher p-value;
* an odds ratio;
* the gene counts for the result table.
The test asks:
Is this phenotype term represented more often in the selected genes than in the rest of the background universe?

It validates the analysis inputs
The function:

.prepare_analysis_inputs()

makes sure that:
* the genes are clean;
* the universe is clean;
* the annotation table is valid;
* selected genes are restricted to the universe;
* genes outside the universe are identified;
* the analysis does not continue with no valid genes.
This fixes an important issue from the earlier prototype, where input genes outside the universe could affect the contingency table.

3. Fisher enrichment
The function:

run_fisher_enrichment()

performs standard enrichment independently for every phenotype term.
For each term, it returns values such as:

selected_with_term
selected_without_term
background_with_term
background_without_term
odds_ratio
p_value
padj

The raw p-values are corrected using:

p.adjust(method = "BH")

This controls the false discovery rate using the Benjamini–Hochberg procedure.
This method provides your basic reference analysis.

4. Parent–child enrichment
The function:

run_parentchild_enrichment()

attempts to reduce the inflation of broad parent terms.
Instead of comparing every term against the entire annotated universe, it compares a term against the genes associated with its direct parent terms.
Conceptually:

Standard Fisher:
specific term versus the entire universe

Parent–child:
specific term versus the relevant parent-term genes

For example, instead of asking:
Is dilated cardiomyopathy enriched compared with all annotated genes?
the method asks something closer to:
Is dilated cardiomyopathy enriched compared with genes already associated with its parent phenotype category?
This can reduce broad and redundant findings.

5. Elim-style enrichment
The function:

run_elim_enrichment()

starts with more specific ontology terms.
When a specific term is significant, its supporting genes are removed from the gene sets of its ancestor terms before those ancestors are tested.
Conceptually:

Specific significant term
        ↓
Remove its supporting genes from ancestors
        ↓
Test broader parent terms

This reduces the possibility that the same group of genes makes every ancestor appear significant.
It should be described as elim-style, because it is inspired by the topGO elim principle but is not necessarily an exact reproduction of every part of the topGO implementation.

6. Optional GO control analysis
The script contains:

run_topgo_control()

This is only for genuine Gene Ontology identifiers beginning with:

GO:

It is not used for HPO or MP.
Its purpose is to allow you to run a GO control analysis when benchmarking the package's ontology-aware methods.
This function is optional and is not needed for the main HPO/MP workflow.

7. It adds readable phenotype names
The function:

add_term_names()

joins phenotype names to the statistical result.
Instead of displaying only:

HP:0001644

the output can show:

HP:0001644 — Dilated cardiomyopathy


8. It prunes redundant terms
The function:

prune_redundant_terms()

checks whether a significant broad term has a significant descendant with strongly overlapping annotated genes.
For example:

Abnormal heart morphology
        ↓
Cardiomyopathy
        ↓
Dilated cardiomyopathy

If the more specific term contains nearly the same supporting genes, the broader ancestor can be removed.
The amount of overlap is controlled using:

prune_overlap = 0.8

This means approximately 80% overlap is required by default.
This is post-analysis pruning. It makes the final results more interpretable, but it does not alter the original ontology or annotation files.

9. It filters unhelpful broad terms
The function:

filter_results_basic()

can remove terms such as:

Mode of inheritance
Autosomal dominant inheritance
Autosomal recessive inheritance
mammalian phenotype

It can also filter terms using:
* maximum background size;
* minimum ontology depth.
This is a simpler filter than ontology pruning and is mainly used to clean report tables and plots.




10. The main user-facing function
The central package function is:

run_pheno_enrichment()

The user supplies:

genes
ontology
method

For example:

result <- run_pheno_enrichment(
  genes = c("LMNA", "MYH7", "TNNT2"),
  ontology = "HPO",
  method = "parentchild"
)

The function then automatically performs:

Clean user genes
        ↓
Load HPO or MP package data
        ↓
Build ontology graph
        ↓
Load annotations
        ↓
Propagate annotations if requested
        ↓
Create or clean the gene universe
        ↓
Remove genes outside the universe
        ↓
Run Fisher, parent–child or elim-style testing
        ↓
Apply FDR correction
        ↓
Add phenotype names
        ↓
Optionally prune redundant terms
        ↓
Return a structured result object


11. HPO and MP selection
The same function supports both ontologies.
HPO

result <- run_pheno_enrichment(
  genes = genes,
  ontology = "HPO"
)

This uses:

human genes
HPO annotations
HPO hierarchy

MP

result <- run_pheno_enrichment(
  genes = genes,
  ontology = "MP"
)

This uses:

human genes
human–mouse one-to-one orthologues
MGI mouse phenotype annotations
MP hierarchy

There are also convenience wrappers:

run_hpo_enrichment()
run_mp_enrichment()


12. What the result object contains
The result is assigned the class:

"pheno_enrichment"

It stores:

ontology
method
annotation source
genes originally supplied
genes included in the analysis
unmapped genes
orthologue mapping summary
universe
universe size
ontology terms
ontology graph
term depths
direct annotations
propagated annotations
raw statistical results
named results
analysis parameters
analysis date and time

This is much better than returning only a data frame because the package retains information needed for:
* summaries;
* plots;
* reproducibility;
* mapping reports;
* method comparisons.

13. Print and summary methods
When the user types:

result

The package prints a concise summary:

Phenotype enrichment
Ontology: HPO
Method: parentchild
Annotation source: HPO
Input genes: 9/10 mapped
Universe size: 4,500
Terms tested: 1,200
Significant FDR<0.05: 15

The function:

summary(result)

provides a more detailed summary and top terms.
This makes the output behave more like a proper R package object.

14. Results tables
The function:

make_results_table()

produces a report-ready table containing information such as:

Ontology
Method
Phenotype
Term ID
Input genes
Gene ratio
Background genes
Odds ratio
P-value
Adjusted p-value
Significant
Overlapping genes

For example:

make_results_table(result, n = 20)

returns the top 20 terms.

15. Visualisations
The script includes:

plot_enrichment()

for a horizontal bar plot, and:

plot_enrichment_bubble()

for a bubble plot.
The plots can show:
* phenotype term names;
* adjusted or raw significance;
* gene ratio;
* number of selected genes;
* the strongest terms.
They are returned as ggplot objects, so users can customise them further.
For example:

p <- plot_enrichment(result)

p + ggplot2::labs(
  title = "Cardiomyopathy HPO enrichment"
)


16. Background-universe comparisons
The function:

compare_backgrounds()

reruns the enrichment using several gene universes.
For example:

backgrounds <- list(
  all_annotated = all_hpo_genes,
  rnaseq_detected = all_detected_genes,
  protein_coding = protein_coding_genes
)

comparison <- compare_backgrounds(
  genes = significant_genes,
  backgrounds = backgrounds,
  ontology = "HPO"
)

This allows you to evaluate how the selected background affects the results, which directly addresses one of your project objectives.

17. Example gene lists
The function:

example_genes()

provides:

example_genes("cardiac")
example_genes("cataract")

These allow users to test the package without preparing their own gene list.

18. Example output support
The script includes:

example_output_table()

to show users the expected result structure.
It also includes:

build_precomputed_example()

which allows me to generate and save a genuine example analysis using your selected ontology and annotation releases.



Biological data handling
* HPO ontology support;
* MP ontology support;
* HPO gene–phenotype annotations;
* MGI mouse phenotype annotations;
* human-to-mouse one-to-one orthologue mapping;
* package-safe processed datasets.


Ontology handling
* OBO parsing;
* parent relationships;
* alternative IDs;
* obsolete-term removal;
* ancestor identification;
* annotation propagation;
* ontology depth calculation.

Statistical analysis
* standard Fisher enrichment;
* Benjamini–Hochberg FDR correction;
* parent–child enrichment;
* elim-style enrichment;
* minimum term-size filtering;
* custom gene universes;
* systematic background comparisons.


Result interpretation
* ontology term names;
* redundant-term pruning;
* broad-term filtering;
* input gene mapping summaries;
* MP orthologue summaries;
* formatted result tables;
* bar plots;
* bubble plots.


Package usability
* one main user function;
* HPO and MP convenience wrappers;
* example gene lists;
* structured result objects;
* print and summary methods;
* package-safe data loading;
* starter unit tests;
* dataset documentation;
* data-preparation scripts.

