#!/usr/bin/bash
 
#SBATCH --job-name=BB_bamindices_merged
#SBATCH --nodes=1
#SBATCH --ntasks=12      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=12:00:00   # modify this to reflect how long to let the job go. 
#SBATCH --output=log_BB_bamindices_merged_%j.txt

module purge
module load python/3.6.5
module load R/3.5.0
module load jdk/18.0.1.1

# Call samtools

samtools="/projects/bwenck@colostate.edu/tools/samtools-1.10/samtools"

$samtools index barcode01_merged_sorted.bam

$samtools index barcode02_merged_sorted.bam

$samtools index barcode03_merged_sorted.bam

$samtools index barcode04_merged_sorted.bam

$samtools index barcode05_merged_sorted.bam

$samtools index barcode06_merged_sorted.bam

$samtools index barcode07_merged_sorted.bam

$samtools index barcode08_merged_sorted.bam

$samtools index barcode09_merged_sorted.bam

$samtools index barcode10_merged_sorted.bam

$samtools index barcode11_merged_sorted.bam

$samtools index barcode12_merged_sorted.bam



