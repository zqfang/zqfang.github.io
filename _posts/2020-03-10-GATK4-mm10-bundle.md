---
title: "GATK4 Mouse Resource Bundle"
description: 'Creating GATK4 mm10 resource bundle'
date: 2020-03-10
categories: ["Make bioinfo uncool again"]
comments: true
---

Collection of GATK4 Mouse Resource Bundle

I found the resource [here](https://github.com/igordot/genomics/blob/master/workflows/gatk-mouse-mm10.md). However, the script is out of date. 
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

### 2. Known dbSNP

Download All in one vcf file from NCBI
```shell
wget ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz \
     -O mouse.dbsnp.vcf.gz
```


Download from the Sanger Mouse Genetics Programme (Sanger MGP)  
{% include alert.html text="Which strain's dbsnp to use?" %}  
Depends on your study.

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

