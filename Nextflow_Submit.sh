#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=32:mem=62gb
#PBS -N Methylation_Submitter_NF
#PBS -j oe

module load nextflow/20.10.0

Project_Dir=/rds/general/user/dthorbur/home/01_Scripts/06_Methylation/07_Nextflow
cd $Project_Dir

echo "Starting: `date`"
nextflow run Methylation.nf -c nextflow.config --profile imperial
echo "All done: `date`"
