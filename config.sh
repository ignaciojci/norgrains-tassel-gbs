#!/bin/bash

#### Basic configuration

RD=/fs/ess/PAS2444/norgrains-tassel-gbs
WD=/fs/ess/PAS2444/25NORGRAINS/production
FASTQ=/fs/ess/PAS2444/24_HistoricalBig6/FastQs
KF=/fs/ess/PAS2444/25NORGRAINS/25NOR01-09_LIBPREPID_03-26-2025.txt
STUDY=25_norgrains_test

#### Advanced configuration (don't usually need to change)
DB=/fs/ess/PAS2444/24_HistoricalBig6/WD/database/HistoricalBig6_032824.db
TASSEL_PL=/fs/ess/PAS2444/TASSEL5/tassel-5-standalone/run_pipeline.pl
RG="/fs/ess/PAS2444/23NORGRAINS/Ensembl_v41_IWGSC_v1.0/Triticum_aestivum.IWGSC.dna.toplevel.fa"
E=PstI-MspI
TAG_LENGTH=75
MIN_Q=0
ext="vcf"
MIN_MEM=100g
MAX_MEM=363g