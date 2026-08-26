# phenoEnrichR

Overview

phenoEnrichR is an R project being developed into an R package for phenotype enrichment analysis. It enables researchers to identify phenotype terms that are significantly overrepresented within a user-supplied human gene set using either the Human Phenotype Ontology (HPO) or the Mammalian Phenotype Ontology (MP).

The project is designed to provide a reproducible workflow for analysing phenotype enrichment and generating publication-ready summary tables and visualisations.

Features
Supports Human Phenotype Ontology (HPO)
Supports Mammalian Phenotype Ontology (MP)
Ontology-aware enrichment analysis
Multiple enrichment methods:
Fisher
Parent-child
Elim
Automatic multiple-testing correction using the Benjamini–Hochberg procedure
Summary tables of enriched phenotype terms
Bar plot and bubble plot visualisations
Analysis Workflow
Input gene list
       │
       ▼
Clean and validate input genes
       │
       ▼
Select ontology
(HPO or MP)
       │
       ▼
Load ontology and phenotype annotations
       │
       ▼
Perform enrichment analysis
       │
       ▼
Adjust p-values (Benjamini–Hochberg)
       │
       ▼
Generate enrichment results
       │
       ▼
Visualise enriched phenotype terms
Example Workflow
Step 1 – Prepare a gene list

Genes can be supplied as:

matrix
data frame
CSV file
TSV file

Example:

genes <- c(
  "LMNA",
  "DSP",
  "MYH7",
  "TTN",
  "PKP2"
)
Step 2 – Run phenotype enrichment
Human Phenotype Ontology
result <- run_pheno_enrichment(
    genes = genes,
    ontology = "HPO",
    method = "parentchild"
)
Mammalian Phenotype Ontology
result <- run_pheno_enrichment(
    genes = genes,
    ontology = "MP",
    method = "parentchild"
)
Step 3 – View the results

Print a summary of the analysis:

print(result)

or

summary(result)
Step 4 – Generate a results table
results_table <- make_results_table(result)

The results table contains:

phenotype name
ontology ID
odds ratio
p-value
adjusted p-value
gene ratio
overlapping genes
Step 5 – Visualise the results
Bar plot
plot_enrichment(result)

Displays the most significantly enriched phenotype terms.

Bubble plot
plot_enrichment_bubble(result)

Displays enriched phenotype terms according to:

statistical significance
gene ratio
number of input genes associated with each phenotype
Supported Enrichment Methods
Method	Description
Fisher	Standard overrepresentation analysis using Fisher's Exact Test.
Parent-child	Accounts for the ontology hierarchy by conditioning each term on its parent terms, reducing enrichment of broad parent phenotypes.
Elim	Tests ontology terms from the most specific to the broadest, reducing redundancy by removing genes contributing to significant child terms before testing ancestor terms.
Supported Ontologies
Human Phenotype Ontology (HPO)

Uses:

Human genes
HPO ontology
Official HPO gene-to-phenotype annotations
Mammalian Phenotype Ontology (MP)

Uses:

Human genes
One-to-one human–mouse orthologues
Mouse Genome Informatics (MGI) phenotype annotations
Mammalian Phenotype Ontology
Output

The enrichment workflow returns a pheno_enrichment object containing:

enrichment statistics
mapped input genes
enriched phenotype terms
adjusted p-values
ontology information
summary statistics
formatted results tables
visualisations

This project is currently being developed into an R package as part of an MSc research project. The core phenotype enrichment workflow has been implemented, including support for HPO and MP enrichment analyses, ontology-aware statistical methods, summary tables, and graphical visualisations.

Overall workflow diagram
┌───────────────────────────────────────────────┐
│ 1. USER INPUT                                 │
│                                               │
│ Human gene list or PanelApp gene panel        │
│ Example: LMNA, DSP, MYH7, TTN, PKP2           │
└──────────────────────┬────────────────────────┘
                       │
                       ▼
┌───────────────────────────────────────────────┐
│ 2. CLEAN AND STANDARDISE GENES                │
│                                               │
│ • Read vector, table, CSV or TSV              │
│ • Select gene column                          │
│ • Remove missing and blank values             │
│ • Convert symbols to uppercase                │
│ • Remove duplicates                           │
└──────────────────────┬────────────────────────┘
                       │
                       ▼
┌───────────────────────────────────────────────┐
│ 3. SELECT ONTOLOGY WORKFLOW                   │
│                                               │
│              HPO             MP               │
│               │               │               │
│ Human genes → HPO      Human genes            │
│ annotations                  │                 │
│                              ▼                 │
│                     Mouse orthologues          │
│                              │                 │
│                              ▼                 │
│                     MGI MP annotations         │
└──────────────────────┬────────────────────────┘
                       │
                       ▼
┌───────────────────────────────────────────────┐
│ 4. LOAD ONTOLOGY STRUCTURE                    │
│                                               │
│ • Term IDs                                    │
│ • Term names                                  │
│ • Parent–child relationships                  │
│ • Ontology depth                              │
│ • Ancestors of each term                      │
└──────────────────────┬────────────────────────┘
                       │
                       ▼
┌───────────────────────────────────────────────┐
│ 5. PROPAGATE ANNOTATIONS                      │
│                                               │
│ A gene annotated to a specific phenotype is   │
│ also assigned to its broader ancestor terms.  │
└──────────────────────┬────────────────────────┘
                       │
                       ▼


┌───────────────────────────────────────────────┐
│ 6. DEFINE THE BACKGROUND UNIVERSE             │
│                                               │
│ Default: approved HGNC protein-coding genes   │
│                                               │
│ Then retain only genes with usable phenotype  │
│ annotations in the selected workflow.         │
└──────────────────────┬────────────────────────┘
                       │
                       ▼
┌───────────────────────────────────────────────┐
│ 7. MAP THE INPUT GENES                        │
│                                               │
│ • Genes used in enrichment                    │
│ • Genes excluded or unmapped                  │
│ • For MP: human-to-mouse mapping summary      │
└──────────────────────┬────────────────────────┘
                       │
                       ▼
┌───────────────────────────────────────────────┐
│ 8. TEST EVERY ELIGIBLE PHENOTYPE TERM         │
│                                               │
│ Choose one method:                            │
│ • Fisher                                      │
│ • Parent–child                                │
│ • Elim                                        │
│                                               │
│ Minimum default term size: 5 genes            │
└──────────────────────┬────────────────────────┘
                       │
                       ▼
┌───────────────────────────────────────────────┐
│ 9. CORRECT FOR MULTIPLE TESTING               │
│                                               │
│ Benjamini–Hochberg adjusted p-values          │
│ Significant by default when FDR ≤ 0.05        │
└──────────────────────┬────────────────────────┘
                       │
                       ▼
┌───────────────────────────────────────────────┐
│ 10. ADD NAMES AND OPTIONALLY PRUNE TERMS      │
│                                               │
│ • Convert IDs to readable phenotype names     │
│ • Optionally remove highly redundant parents  │
└──────────────────────┬────────────────────────┘
                       │
                       ▼
┌───────────────────────────────────────────────┐
│ 11. RETURN RESULTS                            │
│                                               │
│ • Enrichment result object                    │
│ • Summary                                     │
│ • Results table                               │
│ • Enrichment bar plot                         │
│ • Enrichment bubble plot                      │
│ • GEL status plot                             │
│ • GEL mapping plot                            │
└───────────────────────────────────────────────┘


<img width="468" height="642" alt="image" src="https://github.com/user-attachments/assets/0a3dbd6e-cffb-446f-b282-151a93137957" />

