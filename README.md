# PanelAppRex
Aggregating and analysing gene panel data.

## Overview
PanelAppRex is an R project aimed at aggregating and analysing gene panel data. Our repository integrates data that is used genomic testing and the virtual gene panels. This dataset facilitates research and development by providing insights into disease-gene correlations and enhancing variant classification methodologies.

## What you want

PanelAppRex skips the hard work by collecting information about genes and diseases into a simple dataset.

All panels are merged into single tables for your use:
* **Simplified** (Panel ID, Gene)
    - [./data/PanelAppData_combined_minimal.tsv](https://github.com/DylanLawless/PanelAppRex/blob/main/data/PanelAppData_combined_minimal.tsv)
* **Complex** (Panel id, Gene, confidence_level, mode_of_inheritance, name, disease_group, disease_sub_group, status)
    - [./data/PanelAppData_combined_core.tsv](https://github.com/DylanLawless/PanelAppRex/blob/main/data/PanelAppData_combined_core.tsv)

## Example import with R

You can import either the TSV format or Rds format.
Download or clone this repo, then find your preferred format in `/data/`.
The following code is shown in `minimal_example.R`:

```
# TSV format

path_data <- "../data"
path_PanelAppData_genes_combined_core <- paste0(path_data, "/PanelAppData_combined_core")
path_PanelAppData_genes_combined_minimal <- paste0(path_data, "/PanelAppData_combined_minimal")
df_core <- read.table(file= paste0(path_PanelAppData_genes_combined_core, ".tsv"), sep = "\t")
df_minimal <- read.table(file= paste0(path_PanelAppData_genes_combined_minimal, ".tsv"), sep = "\t")

# Rds format

path_data <- "../data"
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/path_PanelAppData_genes_combined_Rds")
df_core <- readRDS(file= path_PanelAppData_genes_combined_Rds)
```

For the full code used see `./src/genomics_england_panels.R`.

## Contents
### Gene panels

**Simplified** (Panel ID, Gene)

```
$ head PanelAppData_combined_minimal.tsv

id      Gene    SYMBOL
1       ABL1    ABL1
1       ACTA2   ACTA2
1       ADAMTSL4        ADAMTSL4
1       ARIH1   ARIH1
1       BGN     BGN
1       COL1A1  COL1A1
1       COL1A2  COL1A2
1       COL3A1  COL3A1
1       COL5A1  COL5A1
```

**Complex** (Panel id, Gene, confidence_level, mode_of_inheritance, name, disease_group, disease_sub_group, status)

```
$ head ./PanelAppData_combined_core.tsv
id      Gene    confidence_level        mode_of_inheritance     name    disease_group   disease_sub_group       status
1       ABL1    3       MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted        Thoracic aortic aneurysm or dissection  Cardiovascular disorders        Connective tissue disorders and aortopathies    public
1       ACTA2   3       MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted        Thoracic aortic aneurysm or dissection  Cardiovascular disorders        Connective tissue disorders and aortopathies    public
1       ADAMTSL4        3       BIALLELIC, autosomal or pseudoautosomal Thoracic aortic aneurysm or dissection  Cardiovascular disorders        Connective tissue disorders and aortopathies    public
1       ARIH1   3       MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted        Thoracic aortic aneurysm or dissection  Cardiovascular disorders        Connective tissue disorders and aortopathies    public
1       BGN     3       X-LINKED: hemizygous mutation in males, monoallelic mutations in females may cause disease (may be less severe, later onset than males) Thoracic aortic aneurysm or dissection Cardiovascular disorders Connective tissue disorders and aortopathies    public
```

### Meta data

```
$ head PanelAppData_combined_meta_names.tsv

panel_id        name
1       Thoracic aortic aneurysm or dissection
3       Stickler syndrome
5       Currarino triad
6       Familial hypercholesterolaemia
...
```

The contents of panels themselves vary in number of genes. 
Some genes are found across many different panels. 
If you are performing WGS analysis and score a variant due to it being a "disease gene", consider that some genes are repeated in many panels and may be unfairly biased.

An example is a panel of 572 PID genes which are well established as consensus in the community.
Another example is a panel with 1675 genes which are associated with unexplained death in infancy and sudden unexplained death in childhood.

## References
The PanelAppRex core model contained 58,592 entries containing annotation fields, including the gene name, disease-gene panel ID, disease-related features, confidence measurements. Data from gnomAD v4 comprised 807,162 individuals, including 730,947 exomes and 76,215 genomes. This dataset provided 786,500,648 single nucleotide variants and 122,583,462 indels, with variant type counts of 9,643,254 synonymous, 16,412,219 missense, 726,924 nonsense, 1,186,588 frameshift and 542,514 canonical splice site variants. ClinVar data were obtained from the variant summary dataset available from the NCBI FTP site, and included 6,845,091 entries, which were processed into 91,319 gene classification groups and a total of 38,983 gene classifications. Data from Ensembl was sourced for validation of identifiers such as gene IDs and Human Genome Organisation Gene Nomenclature Committee (HGNC) symbols. Disease interactions were compared against GE’s PanelApp.

## Recommendations
PanelAppRex is made with care but is produced for research. Clinical applications should refer to accredited clinical sources. 
A widely used source of disease-gene panels is the Genomics England PanelApp, accessible [here](https://panelapp.genomicsengland.co.uk). PanelApp hosts comprehensive gene panels related to genomic tests covered by the NHS, as well as data from historic genomic projects.

## Objectives
The goal of PanelAppRex is to prepare a single dataset that is flexible and instantly usable for human genetic analysis.

## Contributing
Contributions to PanelAppRex are welcome. Please submit a pull request with your updates.

## License
This project is licensed under the MIT License - see the `LICENSE.md` file for details.

## Acknowledgements
Special thanks to Genomics England for providing public access to the PanelApp data.
