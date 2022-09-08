#!/bin/bash


#Call STAR and launch the process for generating the INDEX.

# Job name
#SBATCH --job-name="StarIndex50"

# Runtime and memory
#SBATCH --time=01:59:00

# Partition
#SBATCH --partition=pshort
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --output=./outfile_generated_%j.txt
#SBATCH --error=./errfile_generated_%j.txt

############################################################################

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo -e "\nThe number of threads for this script is: $OMP_NUM_THREADS"


#############################################################################
#### Your shell commands below this line ####

#Load the module for (STAR)
echo -e "Loading the aligner.\n"
module load UHTS/Aligner/STAR/2.7.3a


#Load the GenomeAnalysisToolKit (v4.1)
echo -e "Loading GATK\n"
module load UHTS/Analysis/GenomeAnalysisTK/4.1.0.0

echo -e "Loading samtools\n..."
module load UHTS/Analysis/samtools/1.8


echo -e "Modules loaded"

#############################################################################



########################################
#              PATHS                   #
########################################

#Paths to the data
echo -e "Setting pathways\n"
#Fasta primary assembly.
fastaFiles="/data/users/mbotos/PhDAnalysis/references/ensembl/release102/fasta/Danio_rerio.GRCz11.dna.primary_assembly.fa"
#Annotation gtf file.
#gtfFiles="/data/users/mbotos/benne_andres/references/ensembl/gtfFiles/codingRNAs/Danio_rerio.GRCz11.100.gtf"
gtfFiles="/data/users/mbotos/PhDAnalysis/references/ensembl/release102/gtf/Danio_rerio.GRCz11.102.chr.gtf"

#Create the Index output.
mkdir -p $PWD/primary_assembly_genome_indexed_50 && cd $_

#Set the output of the indexed genome folder.
mkdir -p $PWD/50genomeDirectory
genomeDir="$PWD/50genomeDirectory"



########################################
#         Generate the Index           #
########################################


#Check for SAindexNBases 
#https://www.ncbi.nlm.nih.gov/grc/zebrafish/data
#https://www.wolframalpha.com/input/?i=min%2814%2C+log2%2828.352.955%29+%2F+2+-+1%29

STAR --genomeDir $genomeDir \
--runThreadN $OMP_NUM_THREADS \
--runMode genomeGenerate \
--genomeFastaFiles $fastaFiles \
--sjdbGTFfile $gtfFiles \
--genomeSAindexNbases 10 \
--sjdbOverhang 49


echo -e "Generated the index, now move in the .FAI and .DICT\n"

#GEt out of the first genome and generate the rest.
echo -e "Moving towards last step processes after completing indexing."

#Leave the index genome generation folder.
cd ..

########################################
#             .FA and .DICT            #
########################################

#Reset the name
name=`basename $fastaFiles .fa`
#Set the output path.
out_dict=/data/users/mbotos/PhDAnalysis/references/ensembl/release102/fasta/${name}_50.dict
out_fasta=/data/users/mbotos/PhDAnalysis/references/ensembl/release102/fasta/${name}_50.fai

#Used path "/data/users/mbotos/PhDAnalysis/references/fasta"

################################################
#Create the Dictionary Sequence Reference...   #
################################################

echo -e "Creating Sequence Dictionary\n"
GenomeAnalysisTK CreateSequenceDictionary \
-R=$fastaFiles \
-O=$out_dict


#Completed the dictionary sequence.
echo -e "Completed CreateSequenceDictionary"


#################################################
#Index fasta file.  				#
#################################################
echo -e "Indexing fasta file\n..."		
samtools faidx $fastaFiles 

#As the output just creates the original file in the path specified, take it and change the name the to the one interested.
#mv ${fastaFiles}.fai $out_fasta
#Completed the fasta indexing.
echo -e "Completed fasta indexing"
echo -e "done"
