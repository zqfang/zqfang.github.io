---
title: "Convert Seurat Robj to Scanpy h5ad"
description: The correct way to convert seurat Robj to Scanpy h5ad
date: 2020-04-28
categories: ["Make bioinfo uncool again"]
published: true
comments: true
---

It costed me a lot of time to convert seurat objects to scanpy. It's not a pleasant experience.  
Finally, I solved it.  

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
    adata = r('as.SingleCellExperiment(UpdateSeuratObject(y))')
    adata.write_h5ad(filename=robj.replace("Robj","h5ad"))
```

### 4. other way
seurat -> loom -> scanpy

It's much easier, but I did not test.

1. save to `loom` format fist.
```R
pbmc.loom <- as.loom(pbmc.seurat, filename = "../output/pbmc3k.loom", verbose = FALSE)
pbmc.loom
```
2. read into scanpy
```python
pbmc3k = sc.read_loom("../output/pbmc3k.loom")
```

