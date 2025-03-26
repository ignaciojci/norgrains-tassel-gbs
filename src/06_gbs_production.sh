#!/bin/bash

# ----------------------
# Script for Running TASSEL-GBS Pipeline
# Author: [Your Name]
# Description: This script runs the TASSEL-GBS pipeline with customizable parameters.
# Usage: ./run_tassel_gbs.sh -w /path/to/WD -f /path/to/FastQs -k /path/to/keyfile.txt -s STUDY_NAME -d /path/to/db -r /path/to/reference.fa -e Enzyme -t TAG_LENGTH -q MIN_QUALITY -x OUTPUT_EXT -ms MIN_MEM -mx MAX_MEM
# ----------------------

set -u -e -o pipefail -x

module load bcftools
module load vcftools

# Default values
enzyme="PstI-MspI"
tag_length=75
min_quality=0
output_ext="vcf"
min_mem="100g"
max_mem="300g"

# Parse input arguments
while getopts w:f:k:s:d:r:e:t:q:o:n:x: flag
do
  case "${flag}" in
    w) WD=${OPTARG};;
    f) FASTQ=${OPTARG};;
    k) KF=${OPTARG};;
    s) STUDY=${OPTARG};;
    d) DB=${OPTARG};;
    r) RG=${OPTARG};;
    e) enzyme=${OPTARG};;
    t) tag_length=${OPTARG};;
    q) min_quality=${OPTARG};;
    o) output_ext=${OPTARG};;
    n) min_mem=${OPTARG};;
    x) max_mem=${OPTARG};;
  esac
done

# Check for required parameters
if [[ -z "$WD" || -z "$FASTQ" || -z "$KF" || -z "$STUDY" || -z "$DB" || -z "$RG" ]]; then
  echo "Missing required arguments. Please check your input."
  exit 1
fi

TASSEL_PL=/fs/ess/PAS2444/TASSEL5/tassel-5-standalone/run_pipeline.pl

# Create working directories if not exist
mkdir -p "$WD" "$WD/logs" "$WD/database" "$WD/alignment" "$WD/raw_VCF" "$WD/keyfile"
cd "$WD"
cp "$KF" "$WD/keyfile/"

# Record system information
{
  echo "OS version:"; cat /etc/os-release
  echo -e "\nReference genome & md5 sum:"; md5sum "$RG"
  #echo -e "\nReference genome & md5 sum:"; echo $RG
  echo -e "\nTASSEL version:"; $TASSEL_PL -h 2>&1 | grep "Tassel Version"
  echo -e "\nBCFTools/HTSlib version:"; bcftools --version
  echo -e "\nVCFTools version:"; vcftools --version
} > system_run_info.txt

# Create contig lines for VCF header
if [[ -f "$RG" && ! -f "$RG.fai" ]]; then
  samtools faidx "$RG"
fi
awk 'BEGIN {OFS = ""} {print "##contig=<ID=", $1, ",length=", $2, ">"}' "$RG.fai" > contig_lines.txt

# Output VCF file
OUTFILE="./raw_VCF/${STUDY}_production.${output_ext}"

# Run TASSEL Production
echo "Starting ProductionSNPCallerPluginV2" | tee ./logs/production.log
{
  date
  $TASSEL_PL -Xms${min_mem} -Xmx${max_mem} -fork1 \
    -ProductionSNPCallerPluginV2 -db "$DB" -e "$enzyme" -i "$FASTQ" -k "$KF" \
    -kmerLength "$tag_length" -o "$OUTFILE" -endPlugin -runfork1
} > ./logs/ProductionSNPCallerPluginV2.log

# Capture SLURM job info (if applicable)
if [ ! -z "$SLURM_JOB_ID" ]; then
  scontrol show job=$SLURM_JOB_ID
fi

echo "Pipeline completed. Check logs for details."
