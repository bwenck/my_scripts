#!/usr/bin/bash
 
#SBATCH --job-name=merge_bam_BB
#SBATCH --nodes=1
#SBATCH --ntasks=12      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=12:00:00   # modify this to reflect how long to let the job go. 
#SBATCH --output=log_merge_bam_BB_%j.txt

module purge
module load python/3.6.5
module load R/3.5.0
module load jdk/18.0.1.1

samtools="/projects/bwenck@colostate.edu/tools/samtools-1.10/samtools"

$samtools merge -f ../02_output/aln_reads/bwa/trial/barcode01_merged.bam ../02_output/aln_reads/bwa/trial/FAS71342_pass_barcode01_*.bam

$samtools merge -f ../02_output/aln_reads/bwa/trial/barcode02_merged.bam ../02_output/aln_reads/bwa/trial/FAS71342_pass_barcode02_*.bam

$samtools merge -f ../02_output/aln_reads/bwa/trial/barcode03_merged.bam ../02_output/aln_reads/bwa/trial/FAS71342_pass_barcode03_*.bam

$samtools merge -f ../02_output/aln_reads/bwa/trial/barcode04_merged.bam ../02_output/aln_reads/bwa/trial/FAS71342_pass_barcode04_*.bam

$samtools merge -f ../02_output/aln_reads/bwa/trial/barcode05_merged.bam ../02_output/aln_reads/bwa/trial/FAS71342_pass_barcode05_*.bam

$samtools merge -f ../02_output/aln_reads/bwa/trial/barcode06_merged.bam ../02_output/aln_reads/bwa/trial/FAS71342_pass_barcode06_*.bam

$samtools merge -f ../02_output/aln_reads/bwa/trial/barcode07_merged.bam ../02_output/aln_reads/bwa/trial/FAS71342_pass_barcode07_*.bam

$samtools merge -f ../02_output/aln_reads/bwa/trial/barcode08_merged.bam ../02_output/aln_reads/bwa/trial/FAS71342_pass_barcode08_*.bam

$samtools merge -f ../02_output/aln_reads/bwa/trial/barcode09_merged.bam ../02_output/aln_reads/bwa/trial/FAS71342_pass_barcode09_*.bam

$samtools merge -f ../02_output/aln_reads/bwa/trial/barcode10_merged.bam ../02_output/aln_reads/bwa/trial/FAS71342_pass_barcode10_*.bam

$samtools merge -f ../02_output/aln_reads/bwa/trial/barcode11_merged.bam ../02_output/aln_reads/bwa/trial/FAS71342_pass_barcode11_*.bam

$samtools merge -f ../02_output/aln_reads/bwa/trial/barcode12_merged.bam ../02_output/aln_reads/bwa/trial/FAS71342_pass_barcode12_*.bam
