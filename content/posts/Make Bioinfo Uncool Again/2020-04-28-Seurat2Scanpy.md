---
title: "Convert Seurat to Scanpy h5ad"
description: The correct way to convert seurat to Scanpy h5ad
date: 2020-04-28
categories: ["Make bioinfo uncool again"]
tags: ["scRNA-seq", "Scanpy", "Seurat", "Bioinformatics"]
published: true
comments: true
---

## IMPORTANT UPDATE: 2021-04-15

## SeuratDisk
Please see [SeuratDisk](https://mojaveazure.github.io/seurat-disk/reference/Convert.html) to convert seurat to scanpy.

**Tips**:
1. set default assay to `RNA` before covert to h5ad. 
2. if raw read count need to be imported to anndata, you should only contain counts slot in your seurat object before converstion 


```R
library(SeuratDisk)
# convert factor to character 
i <- sapply(srt@meta.data, is.factor)
srt@meta.data[i] <- lapply(srt@meta.data[i], as.character)
# set default assay
DefaultAssay(srt) <- "RNA"
SaveH5Seurat(srt, filename = "srt.h5seurat", overwrite = TRUE)
Convert("srt.h5seurat", "srt.h5ad", assay="RNA", overwrite = TRUE)
```

## Seurat -> loom -> scanpy

**The best way to convert**: seurat -> loom -> scanpy

It's much easier, since both seurat and scanpy support loom.

1. save to `loom` format.
```R
pbmc.loom <- as.loom(pbmc.seurat.object, filename = "../output/pbmc3k.loom", verbose = FALSE)
pbmc.loom$close_all() # alway close when done 
```
2. read into scanpy
```python
import scanpy as sc
pbmc = sc.read_loom("../output/pbmc.loom", obsm_mapping={"X_umap": ["UMAP_1", "UMAP_2"]})
```

3. open loom in R
```R
immune <- Connect(filename = "../data/mmune_cells.loom", mode = "r")
seurat <- as.Seurat(immune)
immune$close_all() # alway close when done 
```


## rpy2
This is the old way. Very hard to make it work. **Not recommended!**

Convert `Seurat` to `Scanpy` costed me a lot of time to convert seurat objects to scanpy. It's not a pleasant experience.  

### 1. Install `Seurat v3.0.2`, or python kernel will always died!!!
Don't know why latest seurat not work.

### 2. Set the R version for `rpy2`
```python
# user defined R installation
import os
# path to your libR.so, only Seurat v3.0.2 works! 
# create a conda R env for seurat 3.0.2 first
os.environ['R_HOME'] = '/home/fangzq/miniconda/envs/seurat/lib/R' 
# path depends on where you installed Python.
os.environ['R_USER'] = '/home/fangzq/miniconda/lib/python3.7/site-packages/rpy2' 
```
### 3. Now, you'er good to go
```python
import scanpy as sc
import glob
```

Install `anndata2ri` first
```python
import anndata2ri
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
# activate rpy2 env
anndata2ri.activate()
```
```python
robjs = glob.glob("data/*Robj")
```
Convert to h5ad

```python
r('library(Seurat)')
for robj in robjs:
    r(f'x<-load("{robj}")')
    r('y=get(x)')
    r('rm(x)')
    r('DefaultAssay(y) <- "RNA"') # get raw count matrix to save
    # seurat2 object
    # adata = r('as.SingleCellExperiment(UpdateSeuratObject(y))')
    adata = r('as.SingleCellExperiment(y)')
    adata.write_h5ad(filename=robj.replace("Robj","h5ad"))
```



