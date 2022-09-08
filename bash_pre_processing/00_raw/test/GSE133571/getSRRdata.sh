#!/bin/bash

#Download the files from the accesion list with SRA.

# Job name
#SBATCH --job-name="sra"

# Runtime and memory
#SBATCH --time=00:30:00

# Partition
#SBATCH --partition=pall
#SBATCH --mem-per-cpu=6G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-24%12
#SBATCH --output=./outfile_generated_%A_%a_%j.txt
#SBATCH --error=./errfile_generated_%A_%a_%j.txt

############################################################################
 

#Load the modules
echo -e "Loading sra.\n"
module load UHTS/Analysis/sratoolkit/2.10.7
echo -e "Modules loaded"


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo -e "This is array task number $SLURM_ARRAY_TASK_ID"





inputFiles="/data/users/mbotos/PhDdataSets/Datasets/Danio_rerio/RNAseq/GSE133571/SRR_Acc_List.txt"

#Create the array with the links in the file
array=($(<$inputFiles))

########################################
#              PATHS                   #
########################################
#ðownload the data in the following path:
outPath="/data/users/mbotos/PhDdataSets/Datasets/Danio_rerio/RNAseq/GSE133571/"


echo -e "Downloading data...\n Sample: ${array[$SLURM_ARRAY_TASK_ID]}"

########################################
#                RUN                   #
#		 SRA                   #
########################################

fasterq-dump ${array[$SLURM_ARRAY_TASK_ID]} -O $outPath

echo -e "Completed with file ${array[$SLURM_ARRAY_TASK_ID]} as array ID $SLURM_ARRAY_TASK_ID\nDone!"

