#!/bin/bash

#### Basic configuration (all lines)

RD=/fs/ess/PAS2444/norgrains-tassel-gbs # Repository directory
WD=/fs/ess/PAS2444/25NORGRAINS/production
FASTQ=/fs/ess/PAS2444/24_HistoricalBig6/FastQs
KF=/fs/ess/PAS2444/25NORGRAINS/25NOR01-09_LIBPREPID_03-26-2025.txt
STUDY=25_norgrains

#### Basic configuration (IL lines)

RD=/fs/ess/PAS2444/norgrains-tassel-gbs
WD=/fs/ess/PAS2444/25NORGRAINS/test_production_IL
FASTQ=/fs/ess/PAS2444/24_HistoricalBig6/FastQs
KF=/fs/ess/PAS2444/25NORGRAINS/test_IL_GBS_keyfile.txt
STUDY=test_25NOR_IL

#### Advanced configuration (don't usually need to change)
DB=/fs/ess/PAS2444/24_HistoricalBig6/WD/database/HistoricalBig6_032824.db
TASSEL_PL=/fs/ess/PAS2444/TASSEL5/tassel-5-standalone/run_pipeline.pl
RCLONE=/users/PAS1286/jignacio/softwares/rclone-v1.69.1-linux-amd64/rclone
RG="/fs/ess/PAS2444/23NORGRAINS/Ensembl_v41_IWGSC_v1.0/Triticum_aestivum.IWGSC.dna.toplevel.fa"
E=PstI-MspI
TAG_LENGTH=75
MIN_Q=0
ext="vcf"
MIN_MEM=2g
MAX_MEM=363g

#### Scripts
filt_src="$RD/src/filter_markers.sh"