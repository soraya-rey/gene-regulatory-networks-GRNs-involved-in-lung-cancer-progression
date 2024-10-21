# RESYS Project 2023  
**Authors:** Maëlys Delouis, Soraya Rey  
**Institution:** Sorbonne Université, Master 2 Biocomputing and Modeling  

## Overview  
This project was conducted as part of the RESYS unit for the Master 2 program at Sorbonne Université. The goal was to investigate gene regulatory networks (GRNs) involved in lung cancer progression and to identify key regulatory genes across different stages of the disease using multiple network inference algorithms.  

### Main Objectives:  
- Build gene regulatory networks (GRNs) to uncover central genes in lung cancer.  
- Identify both known and novel genes potentially involved in tumor progression.  
- Compare the performance of network inference algorithms: MIIC, GENIE3, and ARACNE.  

---

## Requirements  

### Python  
- Version: 3.10  
- Libraries:  
  - `pandas`, `numpy`, `os`

### R  
- Version: 4.3.2  
- Libraries:  
  - `Seurat`, `GENIE3`, `igraph`, `ggplot2`, `stringr`, `vegan`, `Rgraphviz`, `doParallel`, `doRNG`, `Matrix`, `qgraph`, `minet`, `BiocManager`, `dplyr`, `bnlearn`

### External Tools  
- **MIIC**: [MIIC.curie](https://miic.curie.fr/)  
- **RStudio**  

---

## Project Structure  

- **Report.pdf**: Contains the full project report and bibliography.  
- **/data/**: Contains preprocessed data files.  
  - **Epithelial_cells/**: CSV files with preprocessed data for Normal and Tumoral lung tissue across different cancer stages.  
  - **cell_data_simplified.csv**: Metadata extracted from the original dataset.  

> Original dataset available [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907).

- **/source_code/**:  
  - **Python Scripts:**  
    - `cell_type_separator.py`, `data_process.py`: Scripts to split the original matrix using metadata.  
  - **R Scripts:**  
    - `preprocessing.R`: Seurat pipeline for normalizing data and creating stage-specific datasets.  
    - `network_building.R`: Constructs gene networks using GENIE3 and ARACNE.  
    - `analysis.R`: Generates centrality distributions and identifies key genes in the networks.  
    - `centrality.R`: Compares MIIC centrality and connectivity metrics to classify key regulatory genes.  

- **/results/**:  
  - Contains network matrices (GENIE3, ARACNE, MIIC) and visualizations for each cancer stage, along with centrality distribution plots.  

---

## Running the Project  

1. Update the file paths in the R scripts (`file_names`) to match your local environment.  
2. Open the project in **RStudio** and **VSCode** for execution.  
3. Preprocessing and network-building steps are executed within their respective scripts in the `/source_code/` folder.  

---

## Methods  

### 1. **MIIC (Mutual Information Inference of Causal Networks)**  
- **Description:** MIIC uses mutual information to infer causal gene interactions, capturing both linear and non-linear relationships.  
- **Strengths:** Excels at detecting key interactions and non-linear patterns, particularly useful in biological networks.  
- **Usage:** Applied to discover major regulatory genes at different lung cancer stages.  

### 2. **GENIE3**  
- **Description:** A machine learning algorithm using tree-based models (Random Forest) to infer gene regulatory networks by predicting gene influences.  
- **Strengths:** Handles large datasets and multivariate relationships efficiently.  
- **Usage:** Identified central regulators such as WIF1 and SOX4 in different lung cancer stages.  

### 3. **ARACNE**  
- **Description:** Calculates mutual information to infer direct gene-gene interactions while filtering indirect ones.  
- **Strengths:** Simple and effective for pairwise relationships but sensitive to noise and linear dependencies.  
- **Usage:** Used to find direct interactions, though the results were less reliable than MIIC and GENIE3.  

---

## Data Preprocessing  

1. **Normalization:** Seurat was used to normalize gene expression data, reducing technical variability.  
2. **Filtering:** Genes with low expression or insignificant variance were removed.  
3. **Differential Expression Analysis:** Statistical tests (adjusted t-test) identified genes with significant changes between tumor and normal tissues.  
4. **Gene Selection:** The top 100 genes with the highest log2 fold change (log2FC) were selected for network analysis at each stage.  

---

## Results  

### Key Findings:  

- **Known Genes:** EGFR, KRAS, TP53, BRAF, MYC—well-established in lung cancer progression.  
- **Novel Genes:** Potential new regulators like C9orf24, MS4A8, ITLN2, and RNASE1 were identified and require further validation.  

### Challenges:  
- **Stage-Specific Differences:** Differences in gene networks were observed between cancer stages, particularly in stage 2.  
- **Algorithm Performance:** GENIE3 and MIIC outperformed ARACNE, which struggled with noise and data complexity.  

### Future Directions:  
- Further experimental validation of novel genes and more in-depth exploration of regulatory mechanisms across lung cancer stages.  

