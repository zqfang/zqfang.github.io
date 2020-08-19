---
title: "Deploy snakemake pipeline on HPC"
description: The right way to deploy your snakemake pipeline
date: 2020-08-19
categories: ["Make bioinfo uncool again"]
tags: ["Bioinformatics"]
published: true
comments: true
---

The best part of snakemake is allowed you to run your pipeline on HPC automatically. It save you a lot of time.

## How to run snakemake on HPC
there are two ways to configure

- use `--cluster`: 
  - works on different HPC system, e.g. slurm, SGE. 
  - assign `resource` in `params` directive explicitly.
  - or provide a config file by `--cluster-config` 
- use `--profile`: (recommend way)
  - assign `resource` in `resource` directive explicitly.
  - or provide a profile file.


## Snakemake and slurm

### the `--cluster` way

An example for Stanford Sherlock.

#### 1. define resource in `params` directive

```python
rule eblocks:
    input: "..."
    output: "..."
    params: time = "30:00", mem = "4g" 
    threads: 8
    shell:
       "..."
```

run

```shell
snakemake -s Snakefile --cluster 'sbatch -t {params.time} --mem={params.mem} -c {threads}' -j 10
```

#### 2. define resource in `cluster_config.yaml`

```yaml
# slurm_config.yaml - cluster configuration for Stanford Sherlock
__default__:
    partition: normal
    time_min: "01:00:00" # time limit for each job
    cpus: 1  
    mem: "1g"
    #ntasks-per-node: 14 #Request n cores be allocated per node.
    output: "logs_slurm/{rule}.{wildcards}.out"  ## redirect slurm-JOBID.txt to your directory

strain2trait:
    time_min: "30:00"

eblocks:
    mem: "4g"
    cpus: "{threads}" ## => use `threads` define in rule
```

run 

```shell
snakemake -s Snakefile --cluster-config cluster.yaml \
          --cluster 'sbatch -t {cluster.time_min} --mem={cluster.mem} -c {cluster.cpus} -o {cluster.output} -e {cluster.output}' \
          -j 10
```

#### 3. deploy your pipleline on HPC

make a `submit.sh` script 
```shell
#!/bin/bash
#SBATCH -J snakeflow
#SBATCH --time=120:00:00
#SBATCH --qos long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1  
#SBATCH -p normal
#SBATCH --mem=1g
####SBATCH --mail-type=FAIL
####SBATCH --mail-user=xxx@gmail.com

# activate conda enviroment
source activate base
# run jobs
snakemake -j 666 -s haplomap.smk \
        --configfile config.yaml \
        --cluster "sbatch --time={cluster.time_min} -p {cluster.partition} --mem={cluster.mem} -c {cluster.cpus} " \
        --cluster-config slurm_config.yaml 
```

submit this script to cluster using 

```shell
sbatch submit.sh
```

### the `--profile` way
It's more universal and versatile.

#### 1. create a directory for slurm
```shell
mkdir -p ~/.config/snakemake/slurm
```

#### 2. create `config.yaml` in the slurm directory

**Note**: `resource` directive for clusters only allow integer now. 
```yaml
jobs: 10  # maximun job numbers submitted each time
cluster: "sbatch -p normal -t {resources.time_min} --mem={resources.mem} -c {resources.cpus} -o logs_slurm/{rule}_{wildcards} -e logs_slurm/{rule}_{wildcards} --mail-type=FAIL --mail-user=user@mail.com"
default-resources: [cpus=1, mem=2000, time_min=60]
```

#### 3. assign `resource` if different from default
   
```python
rule eblocks:
    input: "..."
    output: "..."
    resources: mem = 4000 # only allow integers
    shell: "..."
```

#### 4. run 
```shell
snakemake --profile slurm -s haplomap.smk --configfile config.yaml -j 666
```