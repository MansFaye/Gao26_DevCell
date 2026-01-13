#!/usr/bin/env bash

#SBATCH --mail-user=geert.vangeest@bioinformatics.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=24
#SBATCH --output=/data/projects/p682_Single_Cell_RNA_Sequencing_of_MPM/log/infercnv_output_%j.txt
#SBATCH --error=/data/projects/p682_Single_Cell_RNA_Sequencing_of_MPM/log/infercnv_error_%j.txt
#SBATCH --job-name=p682_infercnv
#SBATCH --partition=pall

PROJDIR=/data/projects/p682_Single_Cell_RNA_Sequencing_of_MPM

singularity exec \
--bind $PROJDIR \
$PROJDIR/containers/infercnv.1.3.6-S.simg \
$PROJDIR/cluster/run_inferCNV.R
