#!/bin/bash


#Launch fastp for single end Data.

# Job name
#SBATCH --job-name="fastpS"

# Runtime and memory
#SBATCH --time=00:59:00

# Partition
#SBATCH --partition=pshort
#SBATCH --mem-per-cpu=3G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --array=0-23
#SBATCH --output=./outfile_generated_%A_%a_%j.txt
#SBATCH --error=./errfile_generated_%A_%a_%j.txt

############################################################################
#### Your shell commands below this line ####
############################################################################


#Create output destination.
mkdir -p $PWD/results_fastp_s_nv_2
outputPath="$PWD/results_fastp_s_nv_2"


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#Select thet folder input.
inputFiles="/data/users/mbotos/PhDAnalysis/DataClone/Singlend"


#Get the adaptors list.
adapters="/data/users/mbotos/PhDAnalysis/quality_check/adapters_collection/SingleEnd_adapters.fa"

echo This is array task number $SLURM_ARRAY_TASK_ID

#############################################################################
#module load UHTS/Quality_control/cutadapt/2.3
#module load UHTS/Quality_control/cutadapt/2.5
#module add UHTS/Quality_control/fastp/0.19.5
fastp="/data/users/mbotos/Software/./fastp"


#############################################################################


echo -e "Job starting at: $(date)"

files=($(ls -1 $inputFiles/*.fastq | sort -V))

fq1=${files[$SLURM_ARRAY_TASK_ID]}
#fq2=$(echo $fq1 | sed 's/R1/R2/g')

#If all path included do not use echo, use basename.
#outfq1=$(echo $fq1 | cut -d. -f1)_trim.fastq
#outfq2=$(echo $outfq1 | sed 's/R1/R2/g')

#Name of the reads
outfq1=($(basename $fq1 | cut -d. -f1)_trim.fastq)
#outfq2=($(echo $outfq1 | sed 's/R1/R2/g'))
#outfq2=($(basename $fq2 | cut -d. f1)_trim.fastq)
#out_html=($(basename | cut -d. -f1).html)
out_html="single_end_nv_2.html"
out_json="single_end_nv_2.json"


$fastp --html $outputPath/$out_html \
--json $outputpath/$out_json \
--adapter_fasta $adapters \
-3 \
-5 \
-M 20 \
-q 20 \
-l 25 \
-g \
-x \
-w $SLURM_CPUS_PER_TASK \
-p \
-P 15 \
-i $fq1 \
-o $outputPath/$outfq1

echo -e "Done"
echo -e "Files $fq1 was trimmed by task number $SLURM_ARRAY_TASK_ID and finished at: $(date)"
