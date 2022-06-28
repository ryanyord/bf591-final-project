# Overview
Shiny app to visualize mRNA-Seq expression data profiling Huntington's Disease vs neurologically normal individuals

Link to hosted webapp: [https://ryord.shinyapps.io/mRNA-Seq-Shiny-Visualization/](https://ryord.shinyapps.io/mRNA-Seq-Shiny-Visualization/)

![website](https://i.imgur.com/TeXw1eF.png)

[Video Overview](https://www.youtube.com/watch?v=27cXyJNO2Us&t=7s)

#### Tab 1: Sample Information Exploration
- Summary table to give an overview of each column
- Displays a table of data and histograms of numeric data
#### Tab 2: Counts Matrix Exploration
- Provides a summary table of counts matrix table
- Scatter plot to show user input filtered genes
- Clustered heatmap of counts post-filter
- User selected PCA plot of first 5 PC's
#### Tab 3: Differential Expression
- Plot of user defined x or y axis variables with threshold colored p-adj data points
- Data table of resulting values after filtering by threshold
#### Tab 4: Geneset Enrichment Analysis Exploration
- Normalized enrichment score (NES) plot of up and down-regulated genes
- Table of results with filtering options and download feature
- NES scatterplot with threshold colored p-adj data points


# Files to Utilize App

Files to upload can be downloaded [https://drive.google.com/drive/folders/1eAhvPheWbBvlQoPf1-PHgODRO1qVt_IP?usp=sharing](https://drive.google.com/drive/folders/1eAhvPheWbBvlQoPf1-PHgODRO1qVt_IP?usp=sharing)


|  | File upload for each tab:                                                |
|-------|----------------------------------------------------------------|
| Tab 1 | SraRunTable.csv                                                |
| Tab 2 | TabGSE64810_mlhd_DESeq2_norm_counts_adjust.csv                 |
| Tab 3 | GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.csv |
| Tab 4 | fgsea_results_GSE64810_mlhd_DESeq2_diffexp_DESeq2.csv          |

* Data obtained from [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810)
