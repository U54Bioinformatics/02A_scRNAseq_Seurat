# scRNA-seq analysis using Seurat 
The Seurat package is used for counts data pre-processing, feature selection, clustering and visualization.  

## Pre-requisites and installation 
R (>3.4.0)

```r
install.packages("Seurat")
library(Seurat)
```

## Workflow
1. Load counts data and annotation data (from scRNA-seq pre-processing pipeline)
2. Create seurat object; filter cells 
3. Generate QC data and plots 
4. Filter based on number of features, counts and percent MT
5. Normalize and scale data 
6. Identify variable features 
7. PCA analysis and visualization of components
8. UMAP/TSNE Clustering and visualization 
9. (Optional) Visualize features (markers, cell-type annotations)

### Options/Arguments (different from the default)


### Runtime  
1-2 hours 

## Updates   


## References and documentation 
Reference: [Butler et al. Nature Biotech, 2018](https://doi.org/10.1038/nbt.4096)

Vignette: Satija lab home page (https://satijalab.org/seurat/vignettes.html)
