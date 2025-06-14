#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 24                       # Number of CPU cores
#SBATCH --mem=1024000               # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH --account cattle_gtexs     #project name
#SBATCH -J MASHR                # Name of the job
#SBATCH --output=slurm_%A.out   # STDOUT
#SBATCH --error=slurm_%A.err    # STDERR
#SBATCH -t 24:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #
TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR
#=========================================================================#
#         Your job script                                                 #
#=========================================================================#
export USE_OPENMP=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

strong_file="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/MASHR_Celltype/cell_components/strong_output/strong_pairs.MashR_input.txt.gz"
random_file="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/MASHR_Celltype/cell_components/random_output/nominal_pairs.1000000_subset.MashR_input.txt.gz"

Rscript run_MashR.R ${strong_file} ${random_file} 1 /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/MASHR_Celltype/cell_components/strong_output/output_top_paris

#Rscript run_MashR.R ${strong_file} ${random_file} 0 /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/MASHR_Celltype/Tissue_sharing/all/Final_output1/output_top_paris_across_all

#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
