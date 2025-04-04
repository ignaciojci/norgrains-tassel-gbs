#!/bin/bash

set -u -e -o pipefail -x

config_file="/fs/ess/PAS2444/norgrains-tassel-gbs/config.sh"

source "$config_file"
mkdir -p "$RD/logs" "$RD/tmp"

# Run with
# /fs/ess/PAS2444/norgrains-tassel-gbs/run.sh

# # Detect memory
# if [ -z ${SLURM_MEM_PER_CPU+0} ]; then
#   echo "The environment variable SLURM_MEM_PER_CPU does not exist or is empty." >&2
#   if [ -z ${SLURM_MEM_PER_NODE+0} ]; then
#     echo "The environment variable SLURM_MEM_PER_NODE does not exist or is empty." >&2
#   else
#     echo "The environment variable SLURM_MEM_PER_NODE exists and its value is: ${SLURM_MEM_PER_NODE}." >&2
#     MAX_MEM=${SLURM_MEM_PER_NODE}m
#   fi
# else
#   echo "The environment variable SLURM_MEM_PER_CPU exists and its value is: ${SLURM_MEM_PER_CPU}." >&2
#   mempercpu=${SLURM_MEM_PER_CPU%m}
#   totalmem=$((mempercpu * $SLURM_CPUS_PER_TASK))
#   MAX_MEM="${totalmem}m"
# fi

# Create temporary scripts
TMP_SCRIPT06=$(mktemp --tmpdir="$RD/tmp")
cat > "$TMP_SCRIPT06" <<EOF
#!/bin/bash
module load bcftools/1.16
module load vcftools
"$RD/src/06_gbs_production.sh" -w "$WD" \
  -f "$FASTQ" \
  -k "$KF" \
  -s "$STUDY" \
  -d "$DB" \
  -r "$RG" \
  -e "$E" \
  -t $TAG_LENGTH \
  -q $MIN_Q \
  -o "$ext" \
  -n $MIN_MEM \
  -x $MAX_MEM
EOF

# Run steps
# STEP06_ID=$( \
# sbatch --parsable \
# 	--time=7-00:00:00 \
# 	--nodes=1 --mem=$MAX_MEM --partition=largemem \
# 	--job-name=s06_$STUDY \
# 	--account=PAS2444 \
# 	--output=logs/%x-%A.out \
# 	--error=logs/%x-%A.err \
# 	--wrap="module purge && module load samtools/1.16.1 && bash $TMP_SCRIPT06" )
# echo "Submitted Step 1 (TASSEL-GBS production) with Job ID: $STEP06_ID"

# STEP07_ID=$( \
# 	sbatch --parsable \
# 	#--dependency=afterok:$STEP06_ID \
# 	--time=1-00:00:00 \
# 	--cpus-per-task=4 \
# 	--job-name=s07_$STUDY \
# 	--account=PAS2444 \
# 	--output=logs/%x-%A.out \
# 	--error=logs/%x-%A.err \
# 	--wrap="module purge && module load htslib bcftools/1.16 && bash \"$RD/src/07_reheader_vcf.sh\"" )
# echo "Submitted Step 2 (Reheader VCF) with Job ID: $STEP07_ID"

# STEP08_ID=$( \
# 	sbatch --parsable \
# 	--dependency=afterok:$STEP07_ID \
# 	--time=1-00:00:00 \
# 	--cpus-per-task=2 \
# 	--job-name=s08_$STUDY \
# 	--account=PAS2444 \
# 	--output=logs/%x-%A.out \
# 	--error=logs/%x-%A.err \
# 	--wrap="module purge && module load htslib bcftools/1.16 R/4.4.0-gnu11.2 python/3.9-2022.05 && bash \"$RD/src/08_filter_markers.sh\"" )
# echo "Submitted Step 3 (Filter VCF) with Job ID: $STEP08_ID"

STEP08_ID=$( \
	sbatch --parsable \
	--time=1-00:00:00 \
	--cpus-per-task=2 \
	--job-name=s08_$STUDY \
	--account=PAS2444 \
	--output=logs/%x-%A.out \
	--error=logs/%x-%A.err \
	--wrap="module purge && module load htslib bcftools/1.16 python/3.9-2022.05 R/4.4.0-gnu11.2 && bash \"$RD/src/08_filter_markers.sh\"" )
echo "Submitted Step 3 (Filter VCF) with Job ID: $STEP08_ID"

# Print log files:
echo "Log file paths:"
#echo "logs/s06_${STUDY}-${STEP06_ID}.err"
#echo "logs/s07_${STUDY}-${STEP07_ID}.err"
echo "logs/s08_${STUDY}-${STEP08_ID}.err"

# todo: setup dependencies after each step
### JOBID1=$(sbatch --parsable <other_options> <submission_script>)
### sbatch --dependency=afterok:$JOBID1
# todo: call in pipeline, specifying which step to start

