# Convert Seurat to Scanpy h5ad


## IMPORTANT UPDATE: 2021-04-15

## SeuratDisk
Please see [SeuratDisk](https://mojaveazure.github.io/seurat-disk/reference/Convert.html) to convert seurat to scanpy.

```R
library(Seurat)
library(SeuratDisk)

# step 1: Slim down a Seurat object. So you get raw counts, lognorm counts

seu = DietSeurat(
  srt,
  counts = TRUE, # so, raw counts save to adata.raw.X 
  data = TRUE, # so, log1p counts save to adata.X
  scale.data = FALSE, # set to false, or else will save to adata.X
  features = rownames(srt), # export all genes, not just top highly variable genes
  assays = "RNA",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
  misc = TRUE
)

# step 2: factor to character, or else your factor will be number in adata 
i <- sapply(seu@meta.data, is.factor)
seu@meta.data[i] <- lapply(seu@meta.data[i], as.character)

# step 3: convert 
SaveH5Seurat(seu, filename = "srt.h5seurat", overwrite = TRUE)
Convert("srt.h5seurat", "srt.h5ad", assay="RNA", overwrite = TRUE)
```

```python
# load h5ad
import scanpy as sc
adata = sc.read_h5ad("srt.h5ad")
# save counts to layers
adata.layers['counts'] = adata.raw.X.copy()
adata.layers['log1p'] = data.X.copy()

# you need to scale data for downsream tasks if needed.
```



## Seurat -> loom -> scanpy

You actually neeed additional steps when seurat -> loom -> scanpy

see [here](https://github.com/basilkhuder/Seurat-to-RNA-Velocity)

1. save to `loom` format.
```R
pbmc.loom <- as.loom(seurat_object, filename = "../output/pbmc.loom", verbose = FALSE)
pbmc.loom$close_all() # alway close when done 

write.csv(Cells(seurat_object), file = "cellID_obs.csv", row.names = FALSE)
write.csv(Embeddings(seurat_object, reduction = "umap"), file = "cell_embeddings.csv")
write.csv(seurat_object@meta.data$seurat_clusters, file = "clusters.csv")

```
2. read into scanpy
```python
import scanpy as sc
adata = sc.read_loom("../output/pbmc.loom")

sample_obs = pd.read_csv("cellID_obs.csv")
umap_cord = pd.read_csv("cell_embeddings.csv")
cell_clusters = pd.read_csv("clusters_obs.csv")


# now add metadata to the adata 
# e.g.
adata.obsm['X_umap'] = umap_ordered.values # cellID should be matched first
adata.uns['Cluster_colors'] = ...
...
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




