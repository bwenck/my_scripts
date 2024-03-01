#!/usr/bin/bash

#SBATCH --job-name=netseq_WGS_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=12      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=csu
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=12:00:00   # modify this to reflect how long to let the job go.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bwenck@rams.colostate.edu
#SBATCH --output=log_netseq_WGS_analysis_%j.txt

###This script was created to analyze MinION WGS fastq files for SNPS/indels to verify
###positive strain construction in the archaeal organism, Thermococcus kodakarensis,
###utilizing the Alpine Research Supercomputer cluster at University of Colorado Boulder###

module purge
module load jdk/1.8.0
module load samtools/1.16.1

###Provide path to required programs
fastp="/projects/bwenck@colostate.edu/tools/fastp/fastp"
minimap2="/projects/bwenck@colostate.edu/tools/minimap2/minimap2"
#samtools="/projects/bwenck@colostate.edu/tools/samtools-1.10/samtools"
nanocaller="/projects/bwenck@colostate.edu/tools/NanoCaller/NanoCaller"


file_IDs="W5V_1 W5V_2 W5V_3 W5V_4 W5V_5"

for file_ID in $file_IDs
do
###Call fastp for trimming and quality control

$fastp -i ../01_input/W5V/raw_reads_${file_ID}.fastq.gz -o ../02_output/${file_ID}_trim.fastq -h ../02_output/${file_ID}_report.html -j ../02_output/${file_ID}_report.json --thread ${SLURM_NTASKS} -A -Q -L -G -p

###Call minimap2 to align the filtered reads to the reference genome

$minimap2 -ax map-ont ../01_input/TS559_reference_genome.fasta ../02_output/${file_ID}_trim.fastq > ../02_output/aln_${file_ID}.sam

###Call samtools and convert .sam to .bam, sort, and index

samtools view -Shu ../02_output/aln_${file_ID}.sam > ../02_output/aln_${file_ID}.bam

samtools sort -@ ${SLURM_NTASKS} -o ../02_output/aln_${file_ID}_sorted.bam -O BAM ../02_output/aln_${file_ID}.bam

samtools index -@ ${SLURM_NTASKS} ../02_output/aln_${file_ID}_sorted.bam ../02_output/aln_${file_ID}_sorted.bai

###Activate nanocaller_env and run NanoCaller on sorted.bam files, find SNPS and indels, convert to snps.phased.vcf.gz files

module load anaconda/2022.10
conda activate nanocaller_env

$nanocaller --cpu ${SLURM_NTASKS} --bam ../02_output/aln_${file_ID}_sorted.bam --ref ../01_input/TS559_reference_genome.fasta --preset ont --output ../02_output --prefix ${file_ID}_variants

conda deactivate

done