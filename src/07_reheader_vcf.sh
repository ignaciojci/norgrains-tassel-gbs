#!/bin/bash

#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name=25NORGRAINS_Step7
#SBATCH --account=PAS2444
#SBATCH --output=logs/%x-%A.out
#SBATCH --error=logs/%x-%A.err

# Run with
# sbatch /fs/ess/PAS2444/norgrains-tassel-gbs/src/07_reheader_vcf.sh

set -u -e -o pipefail -x

#### Format the output VCF ####################################################
module load bcftools/1.16
module load vcftools
module load htslib

WD=/fs/ess/PAS2444/25NORGRAINS/production                  ##!!!
STUDY=25_norgrains                                            ##!!!
TASSEL_PL=/fs/ess/PAS2444/TASSEL5/tassel-5-standalone/run_pipeline.pl

source /fs/ess/PAS2444/norgrains-tassel-gbs/config.sh

if [[ -z "${nthreads+x}" ]]; then  # Checks if nthreads is unset
    nthreads="${SLURM_CPUS_PER_TASK:-1}"
fi

# Detect memory
if [ -z ${SLURM_MEM_PER_CPU+0} ]; then
  echo "The environment variable SLURM_MEM_PER_CPU does not exist or is empty." >&2
  if [ -z ${SLURM_MEM_PER_NODE+0} ]; then
    echo "The environment variable SLURM_MEM_PER_NODE does not exist or is empty." >&2
  else
    echo "The environment variable SLURM_MEM_PER_NODE exists and its value is: ${SLURM_MEM_PER_NODE}." >&2
    MAX_MEM=${SLURM_MEM_PER_NODE}m
  fi
else
  echo "The environment variable SLURM_MEM_PER_CPU exists and its value is: ${SLURM_MEM_PER_CPU}." >&2
  mempercpu=${SLURM_MEM_PER_CPU%m}
  totalmem=$((mempercpu * $SLURM_CPUS_PER_TASK))
  MAX_MEM="${totalmem}m"
fi

cd $WD

## Create the name of the output files
OUTFILE=./raw_VCF/$STUDY\_production.vcf
FINAL_FILE=$OUTFILE\.gz

### There's a bunch of stuff in the TASSEL-generated VCF file's header
### lines that can be improved
grep "^##" $OUTFILE > raw_VCF/header.vcf

## If the reference fasta included "chr" before chromosome names, these are
## removed by Tassel, so also need to be stripped out of the contig lines
## in the header for bcftools to work
sed -i 's/ID=chr/ID=/' raw_VCF/header.vcf

## Tassel also capitalizes all characters in chrom names. In the case of wheat
## if the unaligned contigs are given chromosome designator "Un", this must
## be capitalized
sed -i 's/ID=Un/ID=UN/' raw_VCF/header.vcf

## AD and PL FORMAT fields defined incorrectly
## Also, TASSEL's description for the PL field contains equal signs and commas, which causes problems with VCFTools
sed -i 's/AD,Number=./AD,Number=R/g' raw_VCF/header.vcf
sed -i 's/PL,Number=./PL,Number=G/g' raw_VCF/header.vcf
sed -i 's/Description="Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic"/Description="Normalized Phred-scaled likelihoods for genotypes as defined in the VCF specification"/g' raw_VCF/header.vcf

## AF and NS INFO fields do not appear to actually be present in file - delete
grep -E -v "AF,Number=|NS,Number=" raw_VCF/header.vcf > raw_VCF/temp.vcf
mv raw_VCF/temp.vcf raw_VCF/header.vcf

## The "QualityScore" INFO field is not defined in header - add
echo "##INFO=<ID=QualityScore,Number=1,Type=Float,Description=\"Quality Score\">" >> raw_VCF/header.vcf

## Add contig info if present
## Unfortunately Tassel seems to convert any lowercase letters in chrom/contig
## IDs to uppercase, which has caused unexpected headaches
if [[ -f contig_lines.txt ]]; then
cat contig_lines.txt >> raw_VCF/header.vcf
sed -e 's/##contig=<ID=.*,/\U&/' raw_VCF/header.vcf
sed -e 's/##CONTIG/##contig/' raw_VCF/header.vcf
fi

## Update VCF with new header lines
grep -v "^##" $OUTFILE >> raw_VCF/header.vcf
mv raw_VCF/header.vcf $OUTFILE

#### Compress/Index VCF #######################################################
## Compress VCF file to either .vcf.gz or .bcf file
#bcftools view "$OUTFILE" -Oz -o "$FINAL_FILE"
bgzip -c --threads $nthreads "$OUTFILE" > "$FINAL_FILE"
bcftools index -f "$FINAL_FILE"

#### Calculate summary stats ##################################################
echo "Genotype summary" >> ./logs/discovery.log
date >> ./logs/discovery.log
$TASSEL_PL -Xms$MIN_MEM -Xmx$MAX_MEM -vcf "$OUTFILE" -GenotypeSummaryPlugin -endPlugin -export raw_VCF/summary

#### Calculate depth statistics ################################################
vcftools --gzvcf "$OUTFILE" --site-mean-depth --out raw_VCF/raw_VCF

rm "$OUTFILE"
echo "Script finished at:" >> ./logs/discovery.log
date >> ./logs/discovery.log

