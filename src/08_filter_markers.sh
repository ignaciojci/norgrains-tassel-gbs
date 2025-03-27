#!/bin/bash

#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=2
#SBATCH --job-name=25NORGRAINS_Step8
#SBATCH --account=PAS2444
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

# Run with
# sbatch /fs/ess/PAS2444/norgrains-tassel-gbs/src/08_filter_markers.sh

set -u -e -o pipefail -x

module load bcftools/1.16
module load htslib

source /fs/ess/PAS2444/norgrains-tassel-gbs/config.sh
cd "$WD"

if [[ -z "${nthreads+x}" ]]; then  # Checks if nthreads is unset
    nthreads="${SLURM_CPUS_PER_TASK:-1}"
fi

## Create the name of the output files
OUTFILE=./raw_VCF/$STUDY\_production.vcf
FINAL_FILE=$OUTFILE\.gz

vcf_file=$FINAL_FILE

fn=$(basename -- "$vcf_file")

#### Creat output dir
od=${WD}/filt
mkdir -p "$od"

### Remove unknown Chromosomes
echo "Removing unknown chromosome..."
nounchr="$od/04_extracted_dedup_nounchr.vcf.gz"
#zgrep -v ^Un $dedup_file2 | bgzip -c > $nounchr
zgrep -v -i ^un "$vcf_file" | bgzip -c --threads $nthreads > "$nounchr" 
bcftools index -c "$nounchr"

# First, extract all sample names from the VCF
noblanks="$od/04.1_extracted_dedup_nounchr_noblanks.vcf.gz"
bcftools query -l "$nounchr" > "$od/all_samples.txt"

# Filter out unwanted samples (case-insensitive match)
grep -viE "BLANK|NODNA|NO_DNA|NO-DNA|NOTISSUE|NO_TISSUE|NO-TISSUE" "$od/all_samples.txt" > "$od/filtered_samples.txt"

# Create new VCF with only the filtered samples

bcftools view -S "$od/filtered_samples.txt" "$nounchr" -Oz -o "$noblanks" --threads $nthreads
bcftools index -c "$noblanks"

### Filter variants
echo "Filtering..."
#filt_src="/fs/ess/PAS2444/Big6/filt_script/filt_vcf_variants_parallel_CarlosSent_FinalUse.sh"
file_bed=$od/file.bed
bcftools view -h "$noblanks" | grep "contig=" | sed -e 's/^.*ID=//' -e 's/,length=/\t0\t/' -e 's/>//' > "$file_bed"
"$filt_src" "$noblanks" "$od" "$file_bed" 0 $nthreads