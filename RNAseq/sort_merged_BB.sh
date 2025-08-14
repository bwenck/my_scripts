#!/usr/bin/bash
 
#SBATCH --job-name=sort_merged_BB
#SBATCH --nodes=1
#SBATCH --ntasks=12      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=12:00:00   # modify this to reflect how long to let the job go. 
#SBATCH --output=log_sort_merged_BB_%j.txt

module purge
module load python/3.6.5
module load R/3.5.0
module load jdk/18.0.1.1

picard="/projects/bwenck@colostate.edu/tools/picard/build/libs/picard.jar"

java -Djava.io.tmpdir=bla -jar $picard SortSam I=barcode01_merged.bam O=barcode01_merged_sorted.bam SORT_ORDER=coordinate

java -Djava.io.tmpdir=bla -jar $picard SortSam I=barcode02_merged.bam O=barcode02_merged_sorted.bam SORT_ORDER=coordinate

java -Djava.io.tmpdir=bla -jar $picard SortSam I=barcode03_merged.bam O=barcode03_merged_sorted.bam SORT_ORDER=coordinate

java -Djava.io.tmpdir=bla -jar $picard SortSam I=barcode04_merged.bam O=barcode04_merged_sorted.bam SORT_ORDER=coordinate

java -Djava.io.tmpdir=bla -jar $picard SortSam I=barcode05_merged.bam O=barcode05_merged_sorted.bam SORT_ORDER=coordinate

java -Djava.io.tmpdir=bla -jar $picard SortSam I=barcode06_merged.bam O=barcode06_merged_sorted.bam SORT_ORDER=coordinate

java -Djava.io.tmpdir=bla -jar $picard SortSam I=barcode07_merged.bam O=barcode07_merged_sorted.bam SORT_ORDER=coordinate

java -Djava.io.tmpdir=bla -jar $picard SortSam I=barcode08_merged.bam O=barcode08_merged_sorted.bam SORT_ORDER=coordinate

java -Djava.io.tmpdir=bla -jar $picard SortSam I=barcode09_merged.bam O=barcode09_merged_sorted.bam SORT_ORDER=coordinate

java -Djava.io.tmpdir=bla -jar $picard SortSam I=barcode10_merged.bam O=barcode10_merged_sorted.bam SORT_ORDER=coordinate

java -Djava.io.tmpdir=bla -jar $picard SortSam I=barcode11_merged.bam O=barcode11_merged_sorted.bam SORT_ORDER=coordinate

java -Djava.io.tmpdir=bla -jar $picard SortSam I=barcode12_merged.bam O=barcode12_merged_sorted.bam SORT_ORDER=coordinate
