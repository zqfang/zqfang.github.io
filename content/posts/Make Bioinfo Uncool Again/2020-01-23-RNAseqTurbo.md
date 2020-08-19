---
title: "SuperFast RNA-seq"
description: 'salmon-tximport-deseq2'
date: 2020-01-20
categories: ["Make bioinfo uncool again"]
tags: ["RNA-seq","Bioinformatics"]
comments: true
---

salmon-tximport-deseq2
### Step 0: install salmon and download transcriptome cdna from gencode

```shell
conda install salmon
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.transcripts.fa.gz
```

### Step 1. build salmon index
```shell
salmon index -p 8 --gencode -t gencode.v32.transcripts.fa.gz -i salmonIndex_hg38
```

### Step 2: quantification

```shell
salmon quant -i salmonIndex_hg38 -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant
```

### Step 3: merge quantification outputs

use tximport in R
```R
  # R code
  library(tximport)
  library(readr)
  suppressMessages(library('EnsDb.Hsapiens.v86'))

  txdb <- EnsDb.Hsapiens.v86
  k <- keys(txdb, keytype = "GENEID")
  df <- select(txdb, keys = k, keytype = "GENEID", columns = c("TXID","GENEID"))
  tx2gene <- df[, 2:1]  # tx ID, then gene ID

  #tx2gene <- read.table(tx2gene, header= T, sep="\t", stringsAsFactors = F)
  samples <- unlist(strsplit(sample_ids,","))
  salmon.files <- file.path('salmon',samples, "quant.sf")
  names(salmon.files) <- samples
  all(file.exists(salmon.files))
  # get transcript level results
  txi.transcripts <- tximport(salmon.files, type = "salmon", 
                              txOut = TRUE, tx2gene = tx2gene,)
                         #     ignoreTxVersion = TRUE)
  # get gene level results
  txi.salmon <- summarizeToGene(txi.transcripts, tx2gene)

  #save raw counts 
  salmon.counts<- txi.salmon$counts
  salmon.counts<- as.data.frame(salmon.counts)
  write.table(salmon.counts, out_counts, sep="\t", quote=F)
  
  #save gene tpms
  salmon.TPM<- txi.salmon$abundance
  salmon.TPM<- as.data.frame(salmon.TPM)
  write.table(salmon.TPM, out_tpm, sep="\t", quote=F)
  #save transcripts tpms
  salmon.trans.TPM<- txi.transcripts$abundance
  salmon.trans.TPM<- as.data.frame(salmon.trans.TPM)
  write.table(salmon.trans.TPM, outTrans_tpm, sep="\t", quote=F)
   
  save(txi.salmon, file="txi.salmon.RData")

```

### Step 4: Differentially expressed gene analysis

DESeq2 pipeline demo
```R
load("txi.salmon.RData")
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)
dds$condition <- relevel(dds$condition, ref=ctrl)
dds <- DESeq(dds, parallel=TRUE)
res <- results(dds, contrast=c("condition", treat, ctrl))
resOrdered <- res[order(res$padj),]
resOrdered = as.data.frame(resOrdered)
write.table(resOrdered, file="degs.txt", quote=F, sep="\t")
```