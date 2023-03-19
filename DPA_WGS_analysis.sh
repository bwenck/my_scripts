#!/usr/bin/bash
 
#SBATCH --job-name=DPA_WGS_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=12      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=12:00:00   # modify this to reflect how long to let the job go. 
#SBATCH --output=log_DPA_WGS_analysis_%j.txt

module purge
module load python/3.6.5
module load R/3.5.0
module load jdk/1.8.0

#Provide path to required programs
fastp="/projects/bwenck@colostate.edu/tools/fastp/fastp"
bwa="/projects/bwenck@colostate.edu/tools/bwa/bwa"
samtools="/projects/bwenck@colostate.edu/tools/samtools-1.10/samtools"
gatk="/projects/bwenck@colostate.edu/tools/gatk-4.1.7.0/gatk"
picard="/projects/bwenck@colostate.edu/tools/picard/build/libs/picard.jar"


file_IDs="DPA50-1 DPA50-2 DPA50-3 DPA50-4 DPA50-5 DPA50-6 DPA50-7 DPA50-8 DPA50-9 DPA50-10"

for file_ID in $file_IDs
do
#Call fastp for trimming and quality control

$fastp -i ../01_INPUT/DPA_mutants/${file_ID}.1.fastq -I ../01_INPUT/DPA_mutants/${file_ID}.2.fastq -o ../02_OUTPUT/${file_ID}.1_trim.fastq -O ../02_OUTPUT/${file_ID}.2_trim.fastq -h ../02_OUTPUT/${file_ID}_report.html -j ../02_OUTPUT/${file_ID}_report.json --thread ${SLURM_NTASKS} --detect_adapter_for_pe -c -x -g -p

#Call BWA to align the filtered reads to the reference genome

$bwa mem -M -t ${SLURM_NTASKS} ../01_INPUT/DPA_mutants/TS559_v3.fasta ../02_OUTPUT/${file_ID}.1_trim.fastq ../02_OUTPUT/${file_ID}.2_trim.fastq > ../02_OUTPUT/aln_${file_ID}.sam

#Call samtools and convert .sam to .bam

$samtools view -bS ../02_OUTPUT/aln_${file_ID}.sam > ../02_OUTPUT/aln_${file_ID}.bam

#call picard to sort the .bam files

java -Djava.io.tmpdir=bla -jar $picard SortSam I=../02_OUTPUT/aln_${file_ID}.bam O=../02_OUTPUT/aln_${file_ID}_sorted.bam SORT_ORDER=coordinate

$samtools index ../02_OUTPUT/aln_${file_ID}_sorted.bam ../02_OUTPUT/aln_${file_ID}_sorted.bai
done

#Add RG to sorted.bam files

java -jar $picard AddOrReplaceReadGroups I=../02_OUTPUT/aln_DPA50-1_sorted.bam O=../02_OUTPUT/aln_DPA50-1_sorted_RG.bam RGID=HC23T.1 RGLB=lib1 RGPL=ILLUMINA RGPU=HC23TAFX2.1.AGTTCC RGSM=DPA50-1
java -jar $picard AddOrReplaceReadGroups I=../02_OUTPUT/aln_DPA50-2_sorted.bam O=../02_OUTPUT/aln_DPA50-2_sorted_RG.bam RGID=HC23T.1 RGLB=lib2 RGPL=ILLUMINA RGPU=HC23TAFX2.1.ATGTCA RGSM=DPA50-2
java -jar $picard AddOrReplaceReadGroups I=../02_OUTPUT/aln_DPA50-3_sorted.bam O=../02_OUTPUT/aln_DPA50-3_sorted_RG.bam RGID=HC23T.1 RGLB=lib3 RGPL=ILLUMINA RGPU=HC23TAFX2.1.CCGTCC RGSM=DPA50-3
java -jar $picard AddOrReplaceReadGroups I=../02_OUTPUT/aln_DPA50-4_sorted.bam O=../02_OUTPUT/aln_DPA50-4_sorted_RG.bam RGID=HC23T.1 RGLB=lib4 RGPL=ILLUMINA RGPU=HC23TAFX2.1.GTCCGC RGSM=DPA50-4
java -jar $picard AddOrReplaceReadGroups I=../02_OUTPUT/aln_DPA50-5_sorted.bam O=../02_OUTPUT/aln_DPA50-5_sorted_RG.bam RGID=HC23T.1 RGLB=lib5 RGPL=ILLUMINA RGPU=HC23TAFX2.1.GTGAAA RGSM=DPA50-5
java -jar $picard AddOrReplaceReadGroups I=../02_OUTPUT/aln_DPA50-6_sorted.bam O=../02_OUTPUT/aln_DPA50-6_sorted_RG.bam RGID=HC23T.1 RGLB=lib6 RGPL=ILLUMINA RGPU=HC23TAFX2.1.GTGGCC RGSM=DPA50-6
java -jar $picard AddOrReplaceReadGroups I=../02_OUTPUT/aln_DPA50-7_sorted.bam O=../02_OUTPUT/aln_DPA50-7_sorted_RG.bam RGID=HC23T.1 RGLB=lib7 RGPL=ILLUMINA RGPU=HC23TAFX2.1.GTTTCG RGSM=DPA50-7
java -jar $picard AddOrReplaceReadGroups I=../02_OUTPUT/aln_DPA50-8_sorted.bam O=../02_OUTPUT/aln_DPA50-8_sorted_RG.bam RGID=HC23T.1 RGLB=lib8 RGPL=ILLUMINA RGPU=HC23TAFX2.1.CGTACG RGSM=DPA50-8
java -jar $picard AddOrReplaceReadGroups I=../02_OUTPUT/aln_DPA50-9_sorted.bam O=../02_OUTPUT/aln_DPA50-9_sorted_RG.bam RGID=HC23T.1 RGLB=lib9 RGPL=ILLUMINA RGPU=HC23TAFX2.1.GAGTGG RGSM=DPA50-9
java -jar $picard AddOrReplaceReadGroups I=../02_OUTPUT/aln_DPA50-10_sorted.bam O=../02_OUTPUT/aln_DPA50-10_sorted_RG.bam RGID=HC23T.1 RGLB=lib10 RGPL=ILLUMINA RGPU=HC23TAFX2.1.ACTGAT RGSM=DPA50-10

# Run HaplotypeCaller on aln_sorted_RG.bam files, convert to .g.vcf

for file_ID in $file_IDs
do

$samtools index ../02_OUTPUT/aln_${file_ID}_sorted_RG.bam ../02_OUTPUT/aln_${file_ID}_sorted_RG.bai

$gatk --java-options '-Xmx4G' CollectWgsMetrics -I ../02_OUTPUT/aln_${file_ID}_sorted_RG.bam -R ../01_INPUT/DPA_mutants/TS559_v3.fasta -O ../02_OUTPUT/aln_${file_ID}_sorted_RG_wgs_metrics.txt

$gatk --java-options '-Xmx4G' QualityScoreDistribution -I ../02_OUTPUT/aln_${file_ID}_sorted_RG.bam -R ../01_INPUT/DPA_mutants/TS559_v3.fasta -O ../02_OUTPUT/aln_${file_ID}_QSD.txt -CHART ../02_OUTPUT/aln_${file_ID}_QSD.pdf

$gatk --java-options '-Xmx4G' HaplotypeCaller --sample-ploidy=1 \
 -A DepthPerSampleHC \
 -A DepthPerAlleleBySample \
 -A QualByDepth \
 -A MappingQuality \
 -A RMSMappingQuality \
 -A UniqueAltReadCount \
 -A BaseQuality \
 -A BaseQualityRankSumTest \
 -A SampleList \
 -R ../01_INPUT/DPA_mutants/TS559_v3.fasta \
 -I ../02_OUTPUT/aln_${file_ID}_sorted_RG.bam \
 -O ../02_OUTPUT/variants_${file_ID}.g.vcf \
 -ERC GVCF
done
	
#Run .g.vcf files through CombineGVCFs to merge .g.vcf files for GenotypeGVCFs

$gatk --java-options '-Xmx4G' CombineGVCFs -V ../02_OUTPUT/variants_DPA50-1.g.vcf \
 -V ../02_OUTPUT/variants_DPA50-2.g.vcf \
 -V ../02_OUTPUT/variants_DPA50-3.g.vcf \
 -V ../02_OUTPUT/variants_DPA50-4.g.vcf \
 -V ../02_OUTPUT/variants_DPA50-5.g.vcf \
 -V ../02_OUTPUT/variants_DPA50-6.g.vcf \
 -V ../02_OUTPUT/variants_DPA50-7.g.vcf \
 -V ../02_OUTPUT/variants_DPA50-8.g.vcf \
 -V ../02_OUTPUT/variants_DPA50-10.g.vcf \
 -R ../01_INPUT/DPA_mutants/TS559_v3.fasta \
 -O ../02_OUTPUT/DPA_merged_variants.g.vcf

#Run merged.g.vcf files through GenotypeGVCFs to perform joint genotyping

$gatk --java-options '-Xmx4G' GenotypeGVCFs -V ../02_OUTPUT/DPA_merged_variants.g.vcf -R ../01_INPUT/DPA_mutants/TS559_v3.fasta -O ../02_OUTPUT/DPA_joint_variants.vcf

#Split the joint.vcf file into SNPs and indels

java -jar $picard SplitVcfs I=../02_OUTPUT/DPA_joint_variants.vcf SNP_OUTPUT=../02_OUTPUT/DPA_snps.vcf INDEL_OUTPUT=../02_OUTPUT/DPA_indels.vcf STRICT=false

#Filter the variants based on hard filtering conditions

$gatk VariantFiltration \
 -V ../02_OUTPUT/DPA_snps.vcf \
 -filter "QD < 2.0" --filter-name "QD2" \
 -filter "SOR > 3.0" --filter-name "SOR3" \
 -filter "FS > 60.0" --filter-name "FS60" \
 -filter "MQ < 40.0" --filter-name "MQ40" \
 -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
 -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
 -O ../02_OUTPUT/DPA_snps_filtered.vcf
    
$gatk VariantFiltration \
 -V ../02_OUTPUT/DPA_indels.vcf \
 -filter "QD < 2.0" --filter-name "QD2" \
 -filter "FS > 200.0" --filter-name "FS200" \
 -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
 -O ../02_OUTPUT/DPA_indels_filtered.vcf
 
#Make a table of snps and indels for visualization
$gatk VariantsToTable -V ../02_OUTPUT/DPA_snps_filtered.vcf -F POS -F TYPE -F AC -F AF -F AN -F DP -F HOM-REF -F HOM-VAR -GF AD -GF DP -O ../02_OUTPUT/DPA_snps_INFO.table

$gatk VariantsToTable -V ../02_OUTPUT/DPA_indels_filtered.vcf -F POS -F TYPE -F AC -F AF -F AN -F DP -F HOM-REF -F HOM-VAR -GF AD -GF DP -F EVENTLENGTH -O ../02_OUTPUT/DPA_indels_INFO.table