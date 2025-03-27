#!/bin/bash

#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=2
#SBATCH --job-name=25NORGRAINS_Step9
#SBATCH --account=PAS2444
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err 

# Run with
# sbatch /fs/ess/PAS2444/norgrains-tassel-gbs/src/09_split_vcf_by_groups.sh \
#   /fs/ess/PAS2444/25NORGRAINS/production/filt/05_all_regions_filt.vcf.gz \
#   /fs/ess/PAS2444/25NORGRAINS/production/filt/groups.txt \
#   25NORGRAINS 2

set -u -e -o pipefail -x

input_vcf=$1
group_file=$2
prefix=${3:-}
nthreads=${4:-${SLURM_CPUS_PER_TASK:-1}}

# Validate inputs
if [[ ! -f "$input_vcf" ]]; then
    echo "Error: Input VCF file $input_vcf not found" >&2
    exit 1
fi

if [[ ! -f "$group_file" ]]; then
    echo "Error: Group file $group_file not found" >&2
    exit 1
fi

# Load modules
# Try to load bcftools, first specific version then default
if type module &>/dev/null; then
    if module load bcftools/1.16 2>/dev/null; then
        echo "Loaded bcftools/1.16"
    elif module load bcftools 2>/dev/null; then
        echo "Loaded default bcftools"
    fi
fi
module load htslib

source /fs/ess/PAS2444/norgrains-tassel-gbs/config.sh
cd "$WD"

# Create output directory if it doesn't exist
output_dir=$(dirname "$input_vcf")
mkdir -p "$output_dir"

# Filter out unwanted samples (case-insensitive match)
filtered_group_file="${group_file}.noblanks"
grep -viE "BLANK|NODNA|NO_DNA|NO-DNA|NOTISSUE|NO_TISSUE|NO-TISSUE" "$group_file" > "$filtered_group_file"

# Extract unique group names
group_list_file="${group_file}.groups"
awk -F'[ ,\t]' '{print $2}' "$filtered_group_file" | sort -u > "$group_list_file"

while read -r group; do
    # Clean group name for filename
    clean_group=$(echo "$group" | tr -cd '[:alnum:]._-')
    
    # Construct output base name
    if [ -n "$prefix" ] && [[ "$(basename "$input_vcf")" != "${prefix}_"* ]]; then
        base_name="${prefix}_${clean_group}"
    else
        base_name="${clean_group}"
    fi
    
    full_path="${output_dir}/${base_name}"
    
    # Get samples for this group
    sample_list="${full_path}_samples.txt"
    awk -F'[ ,\t]' -v group="$group" '$2 == group {print $1}' "$filtered_group_file" > "$sample_list"
    
    # Skip if no samples found
    if [[ ! -s "$sample_list" ]]; then
        echo "Warning: No samples found for group $group" >&2
        continue
    fi
    
    # Extract these samples from VCF
    if ! bcftools view -S "$sample_list" "$input_vcf" -Oz -o "${full_path}.vcf.gz" --threads "$nthreads"; then
        echo "Error: bcftools view failed for group $group" >&2
        continue
    fi
    
    if ! bcftools index -c "${full_path}.vcf.gz"; then
        echo "Error: bcftools index failed for group $group" >&2
        continue
    fi
    
    echo "Created ${full_path}.vcf.gz for group $group"
    
    # Remove temporary sample list
    rm -f "$sample_list"
done < "$group_list_file"

# Clean up temporary files
#rm -f "$filtered_group_file" "$group_list_file"