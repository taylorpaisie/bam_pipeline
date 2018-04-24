#!/bin/bash
#SBATCH --job-name=bam_pipeline
#SBATCH --output=parallel_bam_pipeline_%j.out
#SBATCH --error=bam_pipeline.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tpaisie@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=80gb
#SBATCH --time=96:00:00
#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --partition=hpg2-compute

pwd; hostname; date

module load intel/2016.0.109 
module load openmpi/1.10.2
module load gcc/5.2.0
module load parallel
module load fastqc
module load trimmomatic
module load bowtie2
module load samtools
module load picard
module load gatk/3.8.0

export _JAVA_OPTIONS="-Xmx8g"

#REF=/ufrc/data/reference/bowtie2/v_cholerae_o1_2010el_1786

#trimmomatic
parallel 'trimmomatic PE {}_1.fastq.gz {}_2.fastq.gz {}_pair_1.fastq.gz {}_unpair_1.fastq.gz {}_pair_2.fastq.gz {}_unpair_2.fastq.gz ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36' ::: $(ls *_*.fastq.gz | rev | cut -c 12- | rev | uniq)

#bowtie2
parallel 'bowtie2 -x /ufrc/data/reference/bowtie2/v_cholerae_o1_2010el_1786 -1 {}_pair_1.fastq.gz -2 {}_pair_2.fastq.gz -I 0 -X 1200 --un-conc {}_unaln_fastq.gz --rg-id "ID_{}" --rg "LB:LB_{}" --rg "PL:ILLUMINA" --rg "SM:SM_{}" -S {}.sam' ::: $(ls *_pair_*.fastq.gz | rev | cut -c 17- | rev | uniq)

# convert sam to bam
parallel 'samtools view -bS {}.sam | samtools sort -o {}_sorted.bam -O BAM' ::: $(ls *.sam | rev | cut -c 5- | rev | uniq)



# Add or replace read groups picard command
# Replace read groups in a BAM file. 
# This tool enables the user to replace all read groups in the INPUT file with a single new read group and assign all reads to this read group in the OUTPUT BAM file.
parallel 'picard AddOrReplaceReadGroups I={}.bam O={}_re.bam RGID=ID_{} RGLB=LB_{} RGPL=ILLUMINA RGPU=ILLUMINA RGSM=SM_{}' ::: $(ls *.bam | rev | cut -c 5- | rev | uniq)

# mark duplicates - Identifies duplicate reads
# This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. 
parallel 'picard MarkDuplicates INPUT={}_re.bam OUTPUT={}_nodups.bam METRICS_FILE=marked_dup_metrics_{}.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true READ_NAME_REGEX="[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*" OPTICAL_DUPLICATE_PIXEL_DISTANCE=100' ::: $(ls *_re.bam | rev | cut -c 8- | rev | uniq)

# index the output bam file from MarkDuplicates
# Must be indexed for the following GATK commands
parallel 'samtools index {}' ::: *_nodups.bam

# my references :)
#REF=/ufrc/salemi/tpaisie/javiana/refseq/NC_020307.fa
#REF=/ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa

# Realigner target creator - Define intervals to target for local realignment
# Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus indel suitable for standard variant discovery approaches
parallel 'GenomeAnalysisTK -T RealignerTargetCreator -R /ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa -I {}_nodups.bam -o forIndelRealigner_{}.intervals' ::: $(ls *_nodups.bam | rev | cut -c 12- | rev | uniq)

# Indel realigner - Perform local realignment of reads around indels
parallel 'GenomeAnalysisTK -T IndelRealigner -R /ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa -I {}_nodups.bam -targetIntervals forIndelRealigner_{}.intervals -o {}_realn.bam' ::: $(ls *_nodups.bam | rev | cut -c 12- | rev | uniq)

# Fix mate information - Verify mate-pair information between mates and fix if needed
# This tool ensures that all mate-pair information is in sync between each read and its mate pair
parallel 'picard FixMateInformation I={}_realn.bam O={}_fix.bam SORT_ORDER=coordinate' ::: $(ls *_realn.bam | rev | cut -c 11- | rev | uniq)

# index final bam file for summary statistics
parallel 'samtools index {}_fix.bam' ::: *_fix.bam


# summary statistics on bam files
# Produces a summary of alignment metrics from a SAM or BAM file
parallel 'picard CollectAlignmentSummaryMetrics R=/ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa I={}_fix.bam O=AlignmentSummaryMetrics_{}.txt ASSUME_SORTED=true MAX_INSERT_SIZE=100000' ::: $(ls *_fix.bam | rev | cut -c 9- | rev | uniq)

# Generate index statistics from a BAM fileThis tool calculates statistics from a BAM index (.bai) file
# The statistics collected include counts of aligned and unaligned reads as well as all records with no start coordinate
parallel 'picard BamIndexStats I={}_fix.bam' ::: $(ls *_fix.bam | rev | cut -c 9- | rev | uniq)

# Assess sequence coverage by a wide array of metrics, partitioned by sample, read group, or library
parallel 'GenomeAnalysisTK -T DepthOfCoverage -R /ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa -I {}_fix.bam -pt sample --outputFormat table' ::: $(ls *_fix.bam | rev | cut -c 9- | rev | uniq)


