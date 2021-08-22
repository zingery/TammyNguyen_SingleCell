#!/bin/bash

module add cellranger/6.0.2

cellranger aggr --id=BM-32-aggr --csv=/project/umw_silvia_corvera/BM/SingleCell_Aug_2021/CellRanger/BM-32-aggr.csv > log_cellranger.log
