#!/bin/bash

module add cellranger/6.0.2

cellranger count --id=BMY-32  --transcriptome=/project/umw_silvia_corvera/software/refdata-gex-GRCh38-2020-A --fastqs=/project/umw_silvia_corvera/BM/Single\ cell_Aug_2021/Fastq --sample=BMY-32  > log_cellranger_bmy.log 
