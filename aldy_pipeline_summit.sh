#!/usr/bin/bash
 
#SBATCH --job-name=mfa_WGS
#SBATCH --nodes=1
#SBATCH --ntasks=12      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=12:00:00   # modify this to reflect how long to let the job go. 
#SBATCH --output=log_mfa_WGS_%j.txt

###This script was designed by breewenck to use in RMACC Summit, University of Colorado, Boulder, CO###

module purge
module load jdk/1.8.0

#Activate conda environment
source /curc/sw/anaconda3/latest
source activate trial_mfa

#Provide path to required programs
fastqc="/projects/geraldy@colostate.edu/tools/FastQC/fastqc"
hisat2="/projects/geraldy@colostate.edu/tools/hisat2/hisat2"
samtools="/projects/geraldy@colostate.edu/tools/samtools-1.10/samtools"
picard="/projects/geraldy@colostate.edu/tools/picard/build/libs/picard.jar"


file_IDs="1 2 3 4 5 6 7 8"

for file_ID in $file_IDs
do

#Run a quality control check on raw gzipped .fastq files before fastp
$fastqc --outdir=/scratch/summit/geraldy@colostate.edu/mfa/02_OUTPUT ${file_ID}.1.fastq.fastqsanger.gz ${file_ID}.2.fastq.fastqsanger.gz 

#Call fastp for trimming and quality control
$fastp -i ../01_INPUT/${file_ID}.1.fastq.fastqsanger.gz -I ../01_INPUT/${file_ID}.2.fastq.fastqsanger.gz -o ../02_OUTPUT/${file_ID}.1_trim.fastq -O ../02_OUTPUT/${file_ID}.2_trim.fastq -h ../02_OUTPUT/${file_ID}_report.html -j ../02_OUTPUT/${file_ID}_report.json --thread ${SLURM_NTASKS} --detect_adapter_for_pe -c -x -g -p -P 1

#Run a quality control check on raw gzipped .fastq files after fastp
$fastqc --outdir=/scratch/summit/geraldy@colostate.edu/mfa/02_OUTPUT ${file_ID}.1_trim.fastq ${file_ID}.2_trim.fastq

#Call HISAT2 to align the filtered reads to the reference genome
$hisat2 --summary-file ../02_OUTPUT/mfa_aln_summary.txt -p ${SLURM_NTASKS} -x ../01_INPUT/TS559 -1 ../02_OUTPUT/${file_ID}.1_trim.fastq -2 ../02_OUTPUT/${file_ID}.2_trim.fastq -S ../02_OUTPUT/aln_${file_ID}.sam

#Call samtools and convert .sam to .bam
$samtools view -bS ../02_OUTPUT/aln_${file_ID}.sam > ../02_OUTPUT/aln_${file_ID}.bam

#call picard to sort the .bam files
java -Djava.io.tmpdir=bla -jar $picard SortSam I=../02_OUTPUT/aln_${file_ID}.bam O=../02_OUTPUT/aln_${file_ID}_sorted.bam SORT_ORDER=coordinate

$samtools index ../02_OUTPUT/aln_${file_ID}_sorted.bam ../02_OUTPUT/aln_${file_ID}_sorted.bai

#Change the bin size as desired
bamCoverage -b ../02_OUTPUT/aln_${file_ID}_sorted.bam -o ../02_OUTPUT/${file_ID}_coverage.bw -bs 100 --normalizeUsing RPKM -e  
done

