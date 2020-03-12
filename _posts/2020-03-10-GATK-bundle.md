---
title: "GATK4 Mouse Resource Bundle"
description: 'Creating GATK4 mm10 resource bundle'
date: 2020-03-10
permalink: /posts/2020/03/blog-post-1/
categories: ["Make bioinfo uncool again"]
comments: true
---

Collection of GATK4 Mouse Resource Bundle

I found the resource [here](https://github.com/igordot/genomics/blob/master/workflows/gatk-mouse-mm10.md). However, the script is out of date. For GATK4, we have

### 1. Genome

Download from NCBI ([mm10](https://www.ncbi.nlm.nih.gov/genome/52)) or Sanger Mouse Genetics Programme
```shell
# NCBI
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz -O GRCm38_68.fa.gz

# or Sanger MGP
wget ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa -O GRCm38_68.fa
```

### 2. Known dbSNP

**Recommend**: Download All in one vcf file from NCBI
```shell
wget ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz \
     -O mouse.dbsnp.vcf.gz
```
or  
Download dbSNP GRCm38 VCF files (each chromosome is in a separate file):
```shell
wget --recursive --no-parent --no-directories \
     --accept vcf.gz \
     ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/
```
then concat
```shell
# generate parameter string containing all VCF files
vcf_file_string=""
for vcf in $(ls -1 vcf_chr_*.vcf) ; do
  vcf_file_string="$vcf_file_string -I $vcf"
done
echo $vcf_file_string

# concatenate VCF files
gatk GatherVcfs -R GRCm38_68.dict $vcf_file_string -O dbsnp.vcf

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
| grep -v "#contig" | grep -v "#source" \
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