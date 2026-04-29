# DEnote
Interactive cell labeling for scRNA-seq data. Built to make a key stage in scRNA-seq analysis accessible to non-bioinformaticians. 

**Disclaimer:** DEnote was initially created while benchmarking / experimenting with [Posit AI](https://posit-ai-beta.share.connect.posit.cloud) interfacing with Claude Sonnet 4.6

## Functionality
- Supports Seurat (.rds) or 10x (.h5) files
- Performs Leiden or Louvain clustering
- Explore 2D/3D UMAP, tSNE, and PCA reductions
- Visualize gene expression across reductions, or violin plots
- De novo marker analysis via FindAllMarkers()
- Score known markers as modules via logistic regression
- Calculate inter-cluster similarity via Bhattacharyya coefficient (BC)
- Manually or automatically label clusters based on visual and statistical criteria
- Export easy to read html reports and .csv of metadata and differential expression results

## User Start Guide
1. Clone github repository
```
git clone https://github.com/Tripfantasy/DEnote.git
```
2. Download dependencies in R (R/Rstudio required) app built using R v4.5.3
```
packages <- c("shiny", "bslib", "Seurat", "plotly", "DT", "ggplot2", 
              "ggrepel", "dplyr", "shinycssloaders", "shinyjs", 
              "rmarkdown", "scales", "RColorBrewer", "tidyr", 
              "tibble", "thematic","bsicons")

install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

lapply(packages, install_if_missing)
```
3. Run App
```
shiny::runApp("path/to/app.R")
```
## Previews

<div align="center">
  <table>
    <tr>
      <td>
        <img width="400" height="300" alt="Screenshot 2026-04-29 at 1 56 08 PM" src="https://github.com/user-attachments/assets/5e50b3da-4f40-4ec7-8600-e2987a5929b9" />
    </td>
    <td>
      <img width="400" height="300" alt="Screenshot 2026-04-29 at 1 48 29 PM" src="https://github.com/user-attachments/assets/4474cb16-6961-4c12-99e8-8994b2da58c9" />
    </td>
  </tr>
</table> 

<div align="center">
  <table>
    <tr>
      <td>
        <img width="400" height="300" alt="Screenshot 2026-04-29 at 1 49 53 PM" src="https://github.com/user-attachments/assets/67776f0c-0c33-4739-93b1-ef60596442cf" />
    </td>
    <td>
      <img width="400" height="300" alt="Screenshot 2026-04-29 at 1 50 24 PM" src="https://github.com/user-attachments/assets/6b8d88aa-5bee-4108-96cd-d6a8d3aa24d9" />
    </td>
  </tr>
</table> 

