#!/bin/bash

#Download the files from the accesion list with SRA.

# Job name
#SBATCH --job-name="sra"

# Runtime and memory
#SBATCH --time=02:30:00

# Partition
#SBATCH --partition=pshort
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=0-2
#SBATCH --output=./outfile_generated_%A_%a_%j.txt
#SBATCH --error=./errfile_generated_%A_%a_%j.txt

############################################################################
 

#Load the modules
echo -e "Loading sra.\n"
module load UHTS/Analysis/sratoolkit/2.10.7
echo -e "Modules loaded"


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo -e "This is array task number $SLURM_ARRAY_TASK_ID"





inputFiles="/data/users/mbotos/PhDdataSets/Datasets/Danio_rerio/RNAseq/Sham/GSE112452/SRR_Acc_List.txt"

#Create the array with the links in the file
array=($(<$inputFiles))

########################################
#              PATHS                   #
########################################
#ðownload the data in the following path:
outPath="/data/users/mbotos/PhDdataSets/Datasets/Danio_rerio/RNAseq/Sham/GSE112452/"


echo -e "Downloading data...\n Sample: ${array[$SLURM_ARRAY_TASK_ID]}"

########################################
#                RUN                   #
#		 SRA                   #
########################################

#fasterq-dump ${array[$SLURM_ARRAY_TASK_ID]} -O $outPath


prefetch ${array[$SLURM_ARRAY_TASK_ID]}

echo -e "Prefetch completed\nStarting dump\n..."
echo -e "test\n $PWD/${array[$SLURM_ARRAY_TASK_ID]}/${array[$SLURM_ARRAY_TASK_ID]}.sra "

fastq-dump --split-e --skip-technical $PWD/${array[$SLURM_ARRAY_TASK_ID]}/${array[$SLURM_ARRAY_TASK_ID]}.sra --outdir $outPath

echo -e "Dumping sra"
echo -e "Completed with file ${array[$SLURM_ARRAY_TASK_ID]} as array ID $SLURM_ARRAY_TASK_ID\nDone!"

