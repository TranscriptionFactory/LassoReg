# #!/bin/bash
# #SBATCH -t 3-00:00
#
# #SBATCH --job-name=runLasso
# # partition (queue) declaration
#
# #SBATCH --mail-type=FAIL
#
# #SBATCH --nodes=1
#
# #SBATCH --ntasks=1
#
# #SBATCH --mem=150g
#
# #SBATCH --cpus-per-task=16
#
# module load gcc/10.2.0
# module load r/4.2.0
#
# # usage: Rscript -d datapath_from_working_dir_including_extension
# Rscript init.R -d data/exampledata.csv
