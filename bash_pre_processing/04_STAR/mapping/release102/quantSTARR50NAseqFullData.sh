#!/bin/bash

#Call STAR quant for each r1-r2 files in the dataset.
#For data to work here and not go over 40 files we select the half of the list so just R1 or R2 and change the mate later with seq if necessary.
#Important note



### --- Mapping with the piRNAclusters our individual datas samples.


# Job name
#SBATCH --job-name="50STARnaseq"

# Runtime and memory
#SBATCH --time=01:59:00

# Partition
#SBATCH --partition=pall
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --array=0-15
##SBATCH --array=0-19%10
#SBATCH --output=./outfile_generated_%A_%a_%j.txt
#SBATCH --error=./errfile_generated_%A_%a_%j.txt
############################################################################


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export file_num=$SLURM_ARRAY_TASK_ID


###########################################################
#	Select the data for mapping			  #
###########################################################

inputFiles="/data/users/mbotos/PhDAnalysis/quality_check/fastp/results_fastp_s_nv_2/results_fastp_s_nv_2/50"

#################################
#   Create output destination   # 
#################################


echo -e "This is array task number $SLURM_ARRAY_TASK_ID"


#############################################################################
#### Your shell commands below this line ####

#Load the module for (STAR)
echo -e "Loading the aligner.\n"
module load UHTS/Aligner/STAR/2.7.3a

echo -e "Modules loaded"

#############################################################################

#######################################
#              PATHS                   #
########################################

#Paths to the data
echo -e "Setting pathways\n"

#Now select indexed Reference
indexedReference="/data/users/mbotos/PhDAnalysis/STAR/indexing/index50/primary_assembly_genome_indexed_50/50genomeDirectory"


#################################################################
#        Create the folder output of all the script.            #
#################################################################


mkdir -p $PWD/sample_$file_num && cd $_
echo -e "Output path created."


#READS Array
#As fastqc takes one by one sample without conisdering pair or single end just proceed with all the files.
files=($(ls -1 $inputFiles/*.fastq)) #| sort -V))

fq1=${files[$SLURM_ARRAY_TASK_ID]}

#fq2=$(echo $fq1 | sed 's/_1/_2/g')

#Name of the reads smaller than 35bp
outfq1=($(basename $fq1 _trim.fastq)) #| cut -d_ -f1,2,5))
alignedOut=${outfq1}_


########################################
#                RUN                   #
#                                      #
#                STAR                  #
#                                      #
#                COUNTS                #
#                                      #
########################################

echo -e "using samples R1 which is $fq1 and \nR2 is $fq2"

STAR --runThreadN $OMP_NUM_THREADS \
--genomeLoad NoSharedMemory \
--runMode alignReads \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMtype BAM Unsorted \
--outFileNamePrefix "$alignedOut" \
--genomeDir $indexedReference \
--readFilesIn $fq1 \
--outFilterMultimapNmax 15 \
--outFilterMismatchNoverLmax 0.06 \
--outFilterMatchNmin 16



#--outFilterScoreMinOverLread 0 \
#--outFilterMatchNminOverLread 0 \
#--limitBAMsortRAM 196000000000 \
#--alignIntronMax 1 \
#--seedSearchStartLmax 20 \
#--alignMatesGapMax 60
#gtf files to access and detect for annotation, get to do them properly.
#--quantMode GeneCounts \
#--sjdbGTFfile $dregff3 \
#--sjdbGTFfile $gtfFiles \
#--sjdbGTFtagExonParentTranscript \
echo -e "Completed Reads Counting"
echo -e "Booooyah"
