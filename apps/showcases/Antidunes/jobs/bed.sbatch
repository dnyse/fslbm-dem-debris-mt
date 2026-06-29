#!/bin/bash
#SBATCH --job-name=antidunes-bed
#SBATCH --partition=singlenode
#SBATCH --constraint=el9
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --time=01:00:00
#SBATCH --output=sbatch-output/output/antidunes_bed_%j.out
#SBATCH --error=sbatch-output/error/antidunes_bed_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dennys.huber@fau.de

# BedGeneration.prm hardcodes blocks <3, 3, 1> = 9 MPI processes
REPO="$HOME/mt-dennys"
"$REPO/apps/showcases/Antidunes/antidunes.sh" --bed --bed-np 9
