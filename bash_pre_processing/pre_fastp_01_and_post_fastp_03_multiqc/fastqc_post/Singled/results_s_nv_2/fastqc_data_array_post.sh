#!/bin/bash


#Call fastqc for each r1-r2 files in the dataset.

# Job name
#SBATCH --job-name="fastqcP"

# Runtime and memory
#SBATCH --time=00:59:00

# Partition
#SBATCH --partition=pshort
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --array=0-23
#SBATCH --output=./outfile_generated_%A_%a_%j.txt
#SBATCH --error=./errfile_generated_%A_%a_%j.txt
############################################################################


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
inputFiles="/data/users/mbotos/PhDAnalysis/quality_check/fastp/results_fastp_s_nv_2/results_fastp_s_nv_2"


#Create output destination.
mkdir -p $PWD/results_s_nv_2
outputPath="$PWD/results_s_nv_2"

echo This is array task number $SLURM_ARRAY_TASK_ID

#############################################################################
#### Your shell commands below this line ####
module load UHTS/Quality_control/fastqc/0.11.7
echo -e "fastqc loaded"

##module load UHTS/Quality_control/cutadapt/2.5
#echo -e "cutadapt loaded"

#############################################################################

#As fastqc takes one by one sample without conisdering pair or single end just proceed with all the files.
#files=($(ls -1 $inputFiles/*.fastq | sort -V))
files=($(ls -1 $inputFiles/*.fastq))


file=${files[$SLURM_ARRAY_TASK_ID]}

##fq2=$(echo $fq1 | sed 's/R1/R2/g')

##If all path included do not use echo, use basename.
##outfq1=$(echo $fq1 | cut -d. -f1)_trim.fastq
##outfq2=$(echo $outfq1 | sed 's/R1/R2/g')



#########################################################################
#Define the name of the file to be renammed by fastqc.

#outfile=($(basename $file | cut -d. -f1))

###For the second does not really makes a difference.
###outfq2=($(echo $outfq1 | sed 's/R1/R2/g'))

#############################################################################

#Now define where the input, what to do with fastqc and where to save it.
echo -e "This is array task number $SLURM_ARRAY_TASK_ID"

echo -e "Working with file:$file"                      

fastqc -t $OMP_NUM_THREADS $file -o $outputPath

echo -e "Started cutadapt with task number $SLURM_ARRAY_TASK_ID\n"
echo -e "Completed with $file\n"
echo -e "File $file was checked, task number is : $SLURM_ARRAY_TASK_ID on $(date)"
