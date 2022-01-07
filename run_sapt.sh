!#/bin/bash
#
#PBS -N mek_sapt0.in
#PBS -l select=1:ncpus=12:mem=64GB
#PBS -q skystd
#PBS -j oe

cd $PBS_O_WORKDIR
module load anaconda3/5.1.0-gcc/8.3.1
source activate mdsapt
python mek_sapt0.py apt_glu_leu_met.yaml mek_res0.csv