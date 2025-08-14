#!/usr/bin/bash
 
#SBATCH --job-name=BB_fC
#SBATCH --nodes=1
#SBATCH --ntasks=12      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=12:00:00   # modify this to reflect how long to let the job go. 
#SBATCH --output=log_BB_fC_%j.txt

module purge
module load python/3.6.5
module load R/3.5.0
module load jdk/18.0.1.1

featureCounts="/projects/bwenck@colostate.edu/tools/subread-2.0.3-Linux-x86_64/bin/featureCounts"

$featureCounts -T ${SLURM_NTASKS} -O -L --verbose -t sncRNA,rRNA,CRISPR,mRNA,tRNA -a ../01_input/TS559_annotations_BRW.gff -G ../01_input/TS559_reference_genome.fasta -o ../02_output/aln_reads/bwa/trial/barcode01_fC_output_2.txt ../02_output/aln_reads/bwa/trial/barcode01_merged_sorted.bam

$featureCounts -T ${SLURM_NTASKS} -O -L --verbose -t sncRNA,rRNA,CRISPR,mRNA,tRNA -a ../01_input/TS559_annotations_BRW.gff -G ../01_input/TS559_reference_genome.fasta -o ../02_output/aln_reads/bwa/trial/barcode02_fC_output_2.txt ../02_output/aln_reads/bwa/trial/barcode02_merged_sorted.bam

$featureCounts -T ${SLURM_NTASKS} -O -L --verbose -t sncRNA,rRNA,CRISPR,mRNA,tRNA -a ../01_input/TS559_annotations_BRW.gff -G ../01_input/TS559_reference_genome.fasta -o ../02_output/aln_reads/bwa/trial/barcode03_fC_output_2.txt ../02_output/aln_reads/bwa/trial/barcode03_merged_sorted.bam

$featureCounts -T ${SLURM_NTASKS} -O -L --verbose -t sncRNA,rRNA,CRISPR,mRNA,tRNA -a ../01_input/TS559_annotations_BRW.gff -G ../01_input/TS559_reference_genome.fasta -o ../02_output/aln_reads/bwa/trial/barcode04_fC_output_2.txt ../02_output/aln_reads/bwa/trial/barcode04_merged_sorted.bam

$featureCounts -T ${SLURM_NTASKS} -O -L --verbose -t sncRNA,rRNA,CRISPR,mRNA,tRNA -a ../01_input/TS559_annotations_BRW.gff -G ../01_input/TS559_reference_genome.fasta -o ../02_output/aln_reads/bwa/trial/barcode05_fC_output_2.txt ../02_output/aln_reads/bwa/trial/barcode05_merged_sorted.bam

$featureCounts -T ${SLURM_NTASKS} -O -L --verbose -t sncRNA,rRNA,CRISPR,mRNA,tRNA -a ../01_input/TS559_annotations_BRW.gff -G ../01_input/TS559_reference_genome.fasta -o ../02_output/aln_reads/bwa/trial/barcode06_fC_output_2.txt ../02_output/aln_reads/bwa/trial/barcode06_merged_sorted.bam

$featureCounts -T ${SLURM_NTASKS} -O -L --verbose -t sncRNA,rRNA,CRISPR,mRNA,tRNA -a ../01_input/TS559_annotations_BRW.gff -G ../01_input/TS559_reference_genome.fasta -o ../02_output/aln_reads/bwa/trial/barcode07_fC_output_2.txt ../02_output/aln_reads/bwa/trial/barcode07_merged_sorted.bam

$featureCounts -T ${SLURM_NTASKS} -O -L --verbose -t sncRNA,rRNA,CRISPR,mRNA,tRNA -a ../01_input/TS559_annotations_BRW.gff -G ../01_input/TS559_reference_genome.fasta -o ../02_output/aln_reads/bwa/trial/barcode08_fC_output_2.txt ../02_output/aln_reads/bwa/trial/barcode08_merged_sorted.bam

$featureCounts -T ${SLURM_NTASKS} -O -L --verbose -t sncRNA,rRNA,CRISPR,mRNA,tRNA -a ../01_input/TS559_annotations_BRW.gff -G ../01_input/TS559_reference_genome.fasta -o ../02_output/aln_reads/bwa/trial/barcode09_fC_output_2.txt ../02_output/aln_reads/bwa/trial/barcode09_merged_sorted.bam

$featureCounts -T ${SLURM_NTASKS} -O -L --verbose -t sncRNA,rRNA,CRISPR,mRNA,tRNA -a ../01_input/TS559_annotations_BRW.gff -G ../01_input/TS559_reference_genome.fasta -o ../02_output/aln_reads/bwa/trial/barcode10_fC_output_2.txt ../02_output/aln_reads/bwa/trial/barcode10_merged_sorted.bam

$featureCounts -T ${SLURM_NTASKS} -O -L --verbose -t sncRNA,rRNA,CRISPR,mRNA,tRNA -a ../01_input/TS559_annotations_BRW.gff -G ../01_input/TS559_reference_genome.fasta -o ../02_output/aln_reads/bwa/trial/barcode11_fC_output_2.txt ../02_output/aln_reads/bwa/trial/barcode11_merged_sorted.bam

$featureCounts -T ${SLURM_NTASKS} -O -L --verbose -t sncRNA,rRNA,CRISPR,mRNA,tRNA -a ../01_input/TS559_annotations_BRW.gff -G ../01_input/TS559_reference_genome.fasta -o ../02_output/aln_reads/bwa/trial/barcode12_fC_output_2.txt ../02_output/aln_reads/bwa/trial/barcode12_merged_sorted.bam