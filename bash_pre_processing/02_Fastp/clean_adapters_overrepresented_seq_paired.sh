#!/bin/bash


#Launch fastp for each r1-r2 files in the dataset, for single end Data.

# Job name
#SBATCH --job-name="fastpP"

# Runtime and memory
#SBATCH --time=00:59:00

# Partition
#SBATCH --partition=pshort
#SBATCH --mem-per-cpu=3G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --array=0-11
#SBATCH --output=./outfile_generated_%A_%a_%j.txt
#SBATCH --error=./errfile_generated_%A_%a_%j.txt

############################################################################
#### Your shell commands below this line ####
############################################################################

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#Select thet folder input.
inputFiles="/data/users/mbotos/PhDAnalysis/DataClone/Pairend"


echo This is array task number $SLURM_ARRAY_TASK_ID

#############################################################################
#module load UHTS/Quality_control/cutadapt/2.3
#module load UHTS/Quality_control/cutadapt/2.5
#module add UHTS/Quality_control/fastp/0.19.5
fastp="/data/users/mbotos/Software/./fastp"


#############################################################################

#Create output destination.
mkdir -p $PWD/results_fastp_p_nv
outputPath="$PWD/results_fastp_p_nv"

#Create the array
files=($(ls -1 $inputFiles/*_1* | sort -V))


#Get the access to each file
fq1=${files[$SLURM_ARRAY_TASK_ID]}
fq2=$(echo $fq1 | sed 's/_1/_2/g')

#If all path included do not use echo, use basename.
#outfq1=$(echo $fq1 | cut -d. -f1)_trim.fastq
#outfq2=$(echo $outfq1 | sed 's/R1/R2/g')

#Name of the output filenames
outfq1=($(basename $fq1 | cut -d. -f1)_trim.fastq)
#outfq2=($(echo $outfq1 | sed 's/R1/R2/g'))
outfq2=($(basename $fq2 | cut -d. -f1)_trim.fastq)
#Name the output html filename
#out_html=($(basename $fq1 | cut -d_ -f1).html)
out_html="pair_end_nv.html"
out_json="pair_end_nv.json"

echo -e "Starting fastp for $fq1\t and $fq2 at: $(date)"
#############################################################################
$fastp -w $SLURM_CPUS_PER_TASK \
--html $outputPath/$out_html \
--json $outputpath/$out_json \
--detect_adapter_for_pe \
-V \
-3 \
-5 \
-M 20 \
-g \
-x \
-p \
-P 20 \
-i $fq1 \
-I $fq2 \
-o $outputPath/$outfq1 \
-O $outputPath/$outfq2


echo -e "Done"
echo -e "Files $fq1 and $fq2 were trimmed by task number $SLURM_ARRAY_TASK_ID and finished at: $(date)"

