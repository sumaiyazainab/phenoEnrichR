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
