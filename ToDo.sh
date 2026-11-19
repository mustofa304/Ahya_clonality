#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=16GB
#SBATCH --mail-type=END
#SBATCH --mail-user=mustofa@usc.edu
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy10_S56_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy11_S57_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy12_S58_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy13_S59_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy14_S60_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy15_S61_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy16_S62_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy17_S63_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy18_S64_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy19_S65_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy20_S66_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy5_S51_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy6_S52_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy7_S53_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy8_S54_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
trim2bRAD_2barcodes_dedup_NoDegen.pl input=Ahy9_S55_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=5000
