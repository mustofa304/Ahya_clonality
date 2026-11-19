#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --mem=16GB
#SBATCH --mail-type=END
#SBATCH --mail-user=mustofa@usc.edu
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

zooxType.pl host= "chr1-12" symgenomes=20-23 > zooxCountsAhya.txt
