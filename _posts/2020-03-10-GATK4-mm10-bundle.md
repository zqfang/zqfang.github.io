---
title: "GATK for genetic analysis of inbred mouse is not good"
description: 
date: 2020-03-10
categories: ["Make bioinfo uncool again"]
published: false
comments: true
---

GATK is design for human genetics, it works poorly on homogeneous inbred mouse.  
So better advoid using it.

One of my colleuge who studies mouse genetics, said, 

> I tried the haplotype caller from GATK. But it seems that the haplotype caller is designed for heterogeneous genome like human than for mice. Therefore, the result coming out of HC is worse than samtools, as I manually inspected a few regions that HC calls didn't make sense.

> In addition, in one of their mouse genomic paper that we reviewed, they even skipped the second recalibration step. We asked them why and they said it was because of the same reason: good for human but not that good for homogeneous inbred mouse.




But we still could collect the resource bundle for mouse.

I found the a workflow [here](https://github.com/igordot/genomics/blob/master/workflows/gatk-mouse-mm10.md). However, the script is out of date. 
Also, see discussion [here](https://www.biostars.org/p/182917/)

For GATK4, we have

### 1. Genome

Download from NCBI ([mm10](https://www.ncbi.nlm.nih.gov/genome/52)) or Sanger Mouse Genetics Programme
```shell
# NCBI
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz -O GRCm38_68.fa.gz

# or Sanger MGP
wget ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa -O GRCm38_68.fa
```

### 2. dbSNP

{% include alert.html text="Which strain's dbsnp to use?" %} 
**Depends on your study design.**


Download All in one vcf file from NCBI
```shell
wget ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz \
     -O mouse.dbsnp.vcf.gz
```

Download from the Sanger Mouse Genetics Programme (Sanger MGP)  

```shell
wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz


```


### 3. Known Indels
For mouse indels, the Sanger Mouse Genetics Programme (Sanger MGP) is probably the best resource.  
Download all MGP indels (5/2015 release):  

```shell
wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz \
-O mgp.v5.indels.vcf.gz
```

Filter for passing variants
```shell
# take header first
zcat mgp.v5.indels.vcf.gz | head -1000 | grep "^#" | cut -f 1-8 \
> mgp.v5.indels.pass.chr.vcf
# keep only passing and append 
zcat mgp.v5.indels.vcf.gz | grep -v "^#" | cut -f 1-8 \
| grep -w "PASS"  >> mgp.v5.indels.pass.chr.vcf
```
Sort VCF (automatically generated index has to be deleted due to a known bug -> No anymore):
```
gatk SortVcf -SD GRCm38_68.dict -I mgp.v5.indels.pass.chr.vcf -O mgp.v5.indels.pass.chr.sort.vcf

# rm .idx
# rm mgp.v5.indels.pass.chr.sort.vcf.idx
```

