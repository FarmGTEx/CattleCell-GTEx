#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 20                       # Number of CPU cores
#SBATCH --mem=102400               # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH --account cattle_gtexs     #project name
#SBATCH -J MASHR                # Name of the job
#SBATCH --output=slurm_%A.out   # STDOUT
#SBATCH --error=slurm_%A.err    # STDERR
#SBATCH -t 10:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #
TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR
#=========================================================================#
#         Your job script                                                 #
#=========================================================================#

###############################################################################################

dir_nominal="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/cell_specific/"
dir_output="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/MASHR_Celltype/Celltype_sharing/random_output/"

### 1. output all top SNP-gene pairs information from permutation results to a .txt file
# (colnames: pheno_id,variant_id,chr,pos)
# file list of permutation results
rm -f -r ${dir_output}/nominal_combined_files.txt
nominal_combined_files=(`find ${dir_nominal} -name "*_LMM.cis_qtl_pairs.*.txt.gz"`)
for l in ${nominal_combined_files[*]}
do
{
    echo ${l} >> ${dir_output}/nominal_combined_files.txt
}
done
./combine_signif_pairs_tjy.py ${dir_output}/nominal_combined_files.txt nominal_pairs -o ${dir_output}
### output file: nominal_pairs.combined_signifpairs.txt.gz
rm -f ${dir_output}/nominal_combined_files.txt


################################################################################################

#tis_names=("Jejunum" "Ileum" "Lung" "MammaryGland" "PituitaryGland" "Placenta" "Testis")
tis_names=($(find "${dir_nominal}" -name "*_LMM.cis_qtl.txt.gz" | \
             sed 's|.*/||' | \
             sed 's/_LMM.cis_qtl.txt.gz//' | \
             sort -u))
             
### 2. extract all nominal pairs from nominal results for each tissue
for ((tis_i=0; tis_i<${#tis_names[@]}; tis_i++))
do
{
    tissue=${tis_names[$tis_i]}
    tissue_name=${tissue%%_*} 
    ct_name=${tissue#*_}
    ct_name=${ct_name#Outliers_}
    nominal_files2=(`ls ${dir_nominal}/${tissue_name}/${ct_name}/${tissue}_LMM.cis_qtl_pairs.*.txt.gz`)
    rm -f ${dir_output}/${tissue}.nominal_files2.txt
    for l in ${nominal_files2[*]}
    do
    {
        echo ${l} >> ${dir_output}/${tissue}.nominal_files2.txt
    }
    done
    # extract_pairs
    python3 ./extract_pairs_tjy.py ${dir_output}/${tissue}.nominal_files2.txt ${dir_output}/nominal_pairs.combined_signifpairs.txt.gz ${tissue}_nominal_pairs -o ${dir_output}
    #> output file: *_nominal_pairs.extracted_pairs.txt.gz
}
done
wait


################################################################################################

#subset_size=10000

### 3. prepare random SNP-gene pairs for MashR
nominal_pairs_files=(`ls ${dir_output}/*_nominal_pairs.extracted_pairs.txt.gz`)
rm -f ${dir_output}/nominal_pairs_files.txt
for l in ${nominal_pairs_files[*]}
do
{
    echo ${l} >> ${dir_output}/nominal_pairs_files.txt
}
done
# MashR format file (z-score)
python3 ./mashr_prepare_input.py ${dir_output}/nominal_pairs_files.txt nominal_pairs.${subset_size}_subset -o ${dir_output} --only_zscore --dropna 
zcat ${dir_output}/nominal_pairs.${subset_size}_subset.MashR_input.txt.gz | sed -e 's/.nominal_pairs_zval//g' | gzip > ${dir_output}/nominal_pairs.${subset_size}_subset.temp.txt.gz
mv ${dir_output}/nominal_pairs.${subset_size}_subset.temp.txt.gz ${dir_output}/nominal_pairs.${subset_size}_subset.MashR_input.txt.gz

################################################################################################

#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID