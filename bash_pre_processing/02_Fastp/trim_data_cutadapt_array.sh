#!/bin/bash


#Call cutadapt and launch the cutadapt process for each r1-r2 files in the dataset, and separate the reads under 15bp and greather than 35bp.

# Job name
#SBATCH --job-name="cutadpt"

# Runtime and memory
#SBATCH --time=01:59:00

# Partition
#SBATCH --partition=pall
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --array=0-19%10
#SBATCH --output=./outfile_generated_%A_%a_%j.txt
#SBATCH --error=./errfile_generated_%A_%a_%j.txt

############################################################################


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
inputFiles="/data/users/mbotos/benne_andres/downloadData/datasets2"



#Create output destination for small reads.
mkdir -p $PWD/results_small_reads_splitted
outputPath1="$PWD/results_small_reads_splitted"

#Create output destination for medium reads.
mkdir -p $PWD/results_medium_reads_splitted
outputPath2="$PWD/results_medium_reads_splitted"


#Create output destination for long reads.
mkdir -p $PWD/results_long_reads_splitted
outputPath3="$PWD/results_long_reads_splitted"


#Create the output destination for the full dataset without separating them by lengths.
mkdir -p $PWD/results_full_reads
outputPath="$PWD/results_full_reads"




echo This is array task number $SLURM_ARRAY_TASK_ID

#############################################################################
#### Your shell commands below this line ####
#module load UHTS/Quality_control/cutadapt/2.3
module load UHTS/Quality_control/cutadapt/2.5

#############################################################################

files=($(ls -1 $inputFiles/*R1* | sort -V))

fq1=${files[$SLURM_ARRAY_TASK_ID]}
fq2=$(echo $fq1 | sed 's/R1/R2/g')

#If all path included do not use echo, use basename.
#outfq1=$(echo $fq1 | cut -d. -f1)_trim.fastq
#outfq2=$(echo $outfq1 | sed 's/R1/R2/g')

#Name of the reads smaller than 35bp
outfq1=($(basename $fq1 | cut -d. -f1)_trim.fastq.gz)
#For the second does not really makes a difference.
outfq2=($(echo $outfq1 | sed 's/R1/R2/g'))


#Short reads out, m<15
m15_outfq1=($(basename $fq1 | cut -d. -f1)_m15.fastq.gz)
m15_outfq2=($(echo $m15_outfq1 | sed 's/R1/R2/g'))


#Longer reads than 35bp. M 35
M35_outfq1=($(basename $fq1 | cut -d. -f1)_M35.fastq.gz)
M35_outfq2=($(echo $M35_outfq1 | sed 's/R1/R2/g'))


#Full reads set.
full_outfq1=($(basename $fq1 | cut -d. -f1)_all_length.fastq.gz)
full_outfq2=($(basename $full_outfq1 | sed 's/R1/R2/g'))



#############################################################################

#Now define where the input, what to do with cutdapt and where to save it.
echo -e "This is array task number $SLURM_ARRAY_TASK_ID"

#illumina rnapcrprimer 5' from blastn : ACTGTAGAACTCTGAACGTGTAGATCT
#ilumina RNAPCRprimer 5': GATCGTCGGACTGTAGAACTCTGAAC ,DpnII: GATCGTCGGACTGTAGAACTCTG
#ilumina index 30 3': TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC
#Illumina Small RNA Adapter 2 (100% over 9bp) :GCACCCGAG obtained from R2 after running cutadapt.
#Illumina Small RNA Adapter 2 (100% over 9bp) :CTCGGGTGC obtained from R1 after running cutadapt.




echo -e "Working with R1:$fq1 and \nand R2:$fq1\n Saving the results in $outputPath"

#cutadapt -j $OMP_NUM_THREADS -a "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC" -a "TGGAATTCTCGGGTGCCAAGG" -a "CTCGGGTGC" -a "TCGTATGCC" -A "GCACCCGAG" -A "GATCGTCGGACTGTAGAACTCTGAAC" -A "GATCGTCGGACTGTAGAACTCTG" -m 15 -M 35 -e 0.1 -O 9 -o $outputPath/$outfq1 -p $outputPath/$outfq2 $fq1 $fq2


#New version
#cutadapt -j $OMP_NUM_THREADS -m 15 -M 35 -e 0.1 -O 9 -o $outputPath/$outfq1 -p $outputPath/$outfq2 $fq1 $fq2

#New version 2.
#cutadapt -j $OMP_NUM_THREADS -b file:adapterListR1.fasta -B file:adapterListR2.fasta --trim-n -m15 -M35 -e 0.1 -q 30 -O 9 --pair-filter=any --too-short-output $outputPath1/$m15_outfq1 --too-short-paired-output $outputPath1/$m15_outfq2 --too-long-output $outputPath3/$M35_outfq1 --too-long-paired-output $outputPath3/$M35_outfq2 -o $outputPath2/$outfq1 -p $outputPath2/$outfq2 $fq1 $fq2



#New version for full length reads without splitting them out.
cutadapt -j $OMP_NUM_THREADS -b file:adapterListR1.fasta -B file:adapterListR2.fasta --trim-n -m15 -e 0.1 -q 30 -O 9 --pair-filter=any -o $outputPath/$full_outfq1 -p $outputPath/$full_outfq2 $fq1 $fq2

echo -e "Started cutadapt with  a=%a A=%A J=%J j=%j task number $SLURM_ARRAY_TASK_ID\n"
echo -e "Completed with $fq1 and $fq2\n"
echo -e "Files $fq1 and $fq2 were trimmed by task number $SLURM_ARRAY_TASK_ID on $(date)"
#echo -e "Completed with $m15_outfq1 and $outfq1 and $M35_outfq1\n"
echo -e "Completetd with $full_outfq1 and $full_outfq2\nBooyah!"
