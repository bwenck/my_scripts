#!/usr/bin/bash
 
#SBATCH --job-name=build_TS559_indices
#SBATCH --nodes=1
#SBATCH --ntasks=12      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=12:00:00   # modify this to reflect how long to let the job go. 
#SBATCH --output=log_build_TS559_indices_%j.txt

###This script was designed by breewenck to use in RMACC Summit, University of Colorado, Boulder, CO###

module purge
module load python/3.6.5
module load R/3.5.0
module load jdk/1.8.0

#Provide path to required programs
hisat2="/projects/geraldy@colostate.edu/tools/hisat2/hisat2"

$hisat2 hisat2-build -p ${SLURM_NTASKS} reference.fa TS559