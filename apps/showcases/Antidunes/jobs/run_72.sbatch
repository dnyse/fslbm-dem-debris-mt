#!/bin/bash
#SBATCH --job-name=antidunes-72
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=72
#SBATCH --time=06:00:00
#SBATCH --output=sbatch-output/output/antidunes_72_%j.out
#SBATCH --error=sbatch-output/error/antidunes_72_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# domainSize <128, 60, 120> / cellsPerBlock <32, 20, 20> = 4x3x6 = 72 blocks
REPO="$HOME/mt-dennys"
"$REPO/apps/showcases/Antidunes/antidunes.sh" --run --np 72
