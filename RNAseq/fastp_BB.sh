#!/usr/bin/bash
 
#SBATCH --job-name=fastp_v2
#SBATCH --nodes=1
#SBATCH --ntasks=12      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=23:00:00   # modify this to reflect how long to let the job go. 
#SBATCH --output=log_fastp_v2_%j.txt
 
fastp="/projects/bwenck@colostate.edu/tools/fastp/fastp"
file_IDs_01="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25"
file_IDs_02="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26"
file_IDs_03="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28"
file_IDs_04="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64"
file_IDs_05="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79"
file_IDs_06="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30"
file_IDs_07="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43"
file_IDs_08="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58"
file_IDs_09="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28"
file_IDs_10="0 1 2 3 4 5 6 7 8 9"
file_IDs_11="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
file_IDs_12="12 13 15 16 28 30 35 42"

# fastp: trim and quality control:
for file_ID_01 in $file_IDs_01
do    
$fastp -i ../01_input/barcode01/FAS71342_pass_barcode01_cce4370c_${file_ID_01}.fastq.gz -o ../02_output/fastp/trial/FAS71342_pass_barcode01_cce4370c_${file_ID_01}_trim.fastq.gz --length_required 15 --length_limit 0 -w ${SLURM_NTASKS} -q 7 -p --failed_out ../02_output/fastp/trial/barcode01_${file_ID_02}_failed.txt
done

for file_ID_02 in $file_IDs_02
do    
$fastp -i ../01_input/barcode02/FAS71342_pass_barcode02_cce4370c_${file_ID_02}.fastq.gz -o ../02_output/fastp/trial/FAS71342_pass_barcode02_cce4370c_${file_ID_02}_trim.fastq.gz --length_required 15 --length_limit 0 -w ${SLURM_NTASKS} -q 7 -p --failed_out ../02_output/fastp/trial/barcode02_${file_ID_02}_failed.txt

done

for file_ID_03 in $file_IDs_03
do    
$fastp -i ../01_input/barcode03/FAS71342_pass_barcode03_cce4370c_${file_ID_03}.fastq.gz -o ../02_output/fastp/trial/FAS71342_pass_barcode03_cce4370c_${file_ID_03}_trim.fastq.gz --length_required 15 --length_limit 0 -w ${SLURM_NTASKS} -q 7 -p --failed_out ../02_output/fastp/trial/barcode03_${file_ID_03}_failed.txt

done
  
for file_ID_04 in $file_IDs_04
do    
$fastp -i ../01_input/barcode04/FAS71342_pass_barcode04_cce4370c_${file_ID_04}.fastq.gz -o ../02_output/fastp/trial/FAS71342_pass_barcode04_cce4370c_${file_ID_04}_trim.fastq.gz --length_required 15 --length_limit 0 -w ${SLURM_NTASKS} -q 7 -p --failed_out ../02_output/fastp/trial/barcode04_${file_ID_04}_failed.txt

done
    
for file_ID_05 in $file_IDs_05
do    
$fastp -i ../01_input/barcode05/FAS71342_pass_barcode05_cce4370c_${file_ID_05}.fastq.gz -o ../02_output/fastp/trial/FAS71342_pass_barcode05_cce4370c_${file_ID_05}_trim.fastq.gz --length_required 15 --length_limit 0 -w ${SLURM_NTASKS} -q 7 -p --failed_out ../02_output/fastp/trial/barcode05_${file_ID_05}_failed.txt

done
   
for file_ID_06 in $file_IDs_06
do    
$fastp -i ../01_input/barcode06/FAS71342_pass_barcode06_cce4370c_${file_ID_06}.fastq.gz -o ../02_output/fastp/trial/FAS71342_pass_barcode06_cce4370c_${file_ID_06}_trim.fastq.gz --length_required 15 --length_limit 0 -w ${SLURM_NTASKS} -q 7 -p --failed_out ../02_output/fastp/trial/barcode06_${file_ID_06}_failed.txt

done
   
for file_ID_07 in $file_IDs_07
do    
$fastp -i ../01_input/barcode07/FAS71342_pass_barcode07_cce4370c_${file_ID_07}.fastq.gz -o ../02_output/fastp/trial/FAS71342_pass_barcode07_cce4370c_${file_ID_07}_trim.fastq.gz --length_required 15 --length_limit 0 -w ${SLURM_NTASKS} -q 7 -p --failed_out ../02_output/fastp/trial/barcode07_${file_ID_07}_failed.txt

done
   
for file_ID_08 in $file_IDs_08
do    
$fastp -i ../01_input/barcode08/FAS71342_pass_barcode08_cce4370c_${file_ID_08}.fastq.gz -o ../02_output/fastp/trial/FAS71342_pass_barcode08_cce4370c_${file_ID_08}_trim.fastq.gz --length_required 15 --length_limit 0 -w ${SLURM_NTASKS} -q 7 -p --failed_out ../02_output/fastp/trial/barcode08_${file_ID_08}_failed.txt

done
    
for file_ID_09 in $file_IDs_09
do    
$fastp -i ../01_input/barcode09/FAS71342_pass_barcode09_cce4370c_${file_ID_09}.fastq.gz -o ../02_output/fastp/trial/FAS71342_pass_barcode09_cce4370c_${file_ID_09}_trim.fastq.gz --length_required 15 --length_limit 0 -w ${SLURM_NTASKS} -q 7 -p --failed_out ../02_output/fastp/trial/barcode09_${file_ID_09}_failed.txt

done
       
for file_ID_10 in $file_IDs_10
do    
$fastp -i ../01_input/barcode10/FAS71342_pass_barcode10_cce4370c_${file_ID_10}.fastq.gz -o ../02_output/fastp/trial/FAS71342_pass_barcode10_cce4370c_${file_ID_10}_trim.fastq.gz --length_required 15 --length_limit 0 -w ${SLURM_NTASKS} -q 7 -p --failed_out ../02_output/fastp/trial/barcode10_${file_ID_10}_failed.txt

done
       
for file_ID_11 in $file_IDs_11
do    
$fastp -i ../01_input/barcode11/FAS71342_pass_barcode11_cce4370c_${file_ID_11}.fastq.gz -o ../02_output/fastp/trial/FAS71342_pass_barcode11_cce4370c_${file_ID_11}_trim.fastq.gz --length_required 15 --length_limit 0 -w ${SLURM_NTASKS} -q 7 -p --failed_out ../02_output/fastp/trial/barcode11_${file_ID_11}_failed.txt

done

for file_ID_12 in $file_IDs_12
do    
$fastp -i ../01_input/barcode12/FAS71342_pass_barcode12_cce4370c_${file_ID_12}.fastq.gz -o ../02_output/fastp/trial/FAS71342_pass_barcode12_cce4370c_${file_ID_12}_trim.fastq.gz --length_required 15 --length_limit 0 -w ${SLURM_NTASKS} -q 7 -p --failed_out ../02_output/fastp/trial/barcode11_${file_ID_12}_failed.txt
           
done
           