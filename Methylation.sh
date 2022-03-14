#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=4:mem=12gb
#PBS -N NF_Methylation_Coordinator
#PBS -j oe

Project_Dir=/rds/general/user/dthorbur/home/01_Scripts/06_Methylation/09_Nextflow3/FAK69425
cd $Project_Dir

module load nextflow/20.10.0

echo "Starting:`date`"

nextflow run Methylation.nf \
	-c Methylation.config \
	--profile imperial \
	--RunID "FAK69425" \
	--Input "/rds/general/user/dthorbur/home/01_Scripts/06_Methylation/08_Nextflow2/FAK69425/00_Tarball/20191003_1351_MN22778_FAK69425_695686d4.tar.gz" \
	--RefGen "/rds/general/user/dthorbur/home/01_Scripts/06_Methylation/07_Nextflow/00_Resources/VectorBase-54_AgambiaePEST_Genome.fasta" \
	--Container "/rds/general/user/dthorbur/home/01_Scripts/06_Methylation/07_Nextflow/00_Resources/ONT_Guppy_GPU.sif" \
	--Skip_Decompress

echo "Finished:`date`"
