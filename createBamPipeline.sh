#!/usr/bin/env bash
#SBATCH --job-name=bam_pipeline_CH16
#SBATCH --mail-type=END,ABORT
#SBATCH --mail-user=tpaisie@ufl.edu
#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --output=out_createBamPipeline_test_log_%j
#SBATCH --nodes=1
#SBATCH --ntasks=2
##SBATCH --constraint=c6145
#SBATCH --mem=10gb
#SBATCH --time=100:00:00

module load python
module load sickle
module load fastqc
module load bowtie2
module load samtools
module load picard/2.5.0
module load gatk/3.5.0


# converts fastq files to fastq sanger format
python convertFastqToFastqSanger.py CH16

sickle pe -f CH16_1s.fastq -r CH16_2s.fastq -t sanger -o CH16_1_trimmed.fastq -p CH16_2_trimmed.fastq -s CH16_singles.fastq -q 30 -l 30

fastqc -o /ufrc/salemi/tpaisie/cholera/fastqc_files -f fastq CH16_1_trimmed.fastq CH16_2_trimmed.fastq

bowtie2 -x /ufrc/data/reference/bowtie2/v_cholerae_o1_2010el_1786 -1 CH16_1_trimmed.fastq -2 CH16_2_trimmed.fastq -I 0 -X 1200 -S CH16.sam --un-conc unaligned_CH16.sam --rg-id "ID_CH16" --rg "LB:LB_CH16" --rg "PL:ILLUMINA" --rg "SM:SM_CH16"

samtools view -bS CH16.sam | samtools sort -o sorted_CH16.bam -O BAM


# running picard on bam file
export _JAVA_OPTIONS="-Xms1g -Xmx8g"

picard AddOrReplaceReadGroups I=sorted_CH16.bam O=CH16_aorrg.bam RGID=RD_CH16 RGLB=LB_CH16 RGPL=illumina RGPU=PU_CH16 RGSM=SM_CH16 SORT_ORDER=coordinate

picard MarkDuplicates INPUT=CH16_aorrg.bam OUTPUT=dups_CH16.bam METRICS_FILE=marked_dup_metrics_CH16.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true READ_NAME_REGEX="[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" OPTICAL_DUPLICATE_PIXEL_DISTANCE=100

samtools index dups_CH16.bam


GenomeAnalysisTK -T RealignerTargetCreator -R /ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa -I dups_CH16.bam -o forIndelRealigner.intervals

GenomeAnalysisTK -T IndelRealigner -R /ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa -I dups_CH16.bam -targetIntervals forIndelRealigner.intervals -o realign_CH16.bam

picard FixMateInformation I=realign_CH16.bam O=fix_CH16.bam SORT_ORDER=coordinate

samtools index fix_CH16.bam 

# summary statistics

picard CollectAlignmentSummaryMetrics R=/ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa I=fix_CH16.bam O=AlignmentSummaryMetrics_CH16.txt ASSUME_SORTED=true MAX_INSERT_SIZE=100000

picard BamIndexStats I=fix_CH16.bam

GenomeAnalysisTK -T DepthOfCoverage -R /ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa -I fix_CH16.bam -pt sample --outputFormat table

