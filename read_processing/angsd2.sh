#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=16GB
#SBATCH --mail-type=END
#SBATCH --mail-user=mustofa@usc.edu
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

module load htslib
module load zlib
module load bzip2
module load curl

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -hwe_pval 1e-5 -skipTriallelic 1 -minInd 189 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doDepth 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b bams2 -GL 1 $FILTERS $TODO -P 1 -out ahya2

