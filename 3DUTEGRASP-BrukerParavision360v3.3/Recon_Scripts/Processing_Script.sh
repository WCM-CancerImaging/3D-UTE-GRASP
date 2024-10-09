#! /bin/bash -l

#SBATCH --partition=scu-gpu   # cluster-specific: GPU nodes

#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --job-name=SCU_Processing_Script
#SBATCH --time=48:00:00   # HH/MM/SS
#SBATCH --mem=30G   # memory requested, units available: K,M,G,T
#SBATCH --output=oRecon_GRASP_PV360.out

#SBATCH --gres=gpu:1 # request GPU

#spack load cuda@10.2.89%gcc@4.8.5

cd /athena/kimlab/scratch/jzh4009/Recon_3DUTEGRASP_PV360/Recon_Scripts
cat Processing_Script.m | srun /home/software/apps/matlab/R2020a/bin/matlab -nodisplay
