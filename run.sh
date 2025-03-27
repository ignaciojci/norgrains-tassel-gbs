#!/bin/bash

#SBATCH --time=7-00:00:00
#SBATCH --nodes=1 --mem=363g --partition=largemem
#SBATCH --job-name=25NORGRAINS
#SBATCH --account=PAS2444
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

set -u -e -o pipefail -x

module load bcftools/1.16
module load vcftools

# Run with
# sbatch /fs/ess/PAS2444/norgrains-tassel-gbs/run.sh

config_file="/fs/ess/PAS2444/norgrains-tassel-gbs/config.sh"

source "$config_file"

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

# todo: setup dependencies after each step
### JOBID1=$(sbatch --parsable <other_options> <submission_script>)
### sbatch --dependency=afterok:$JOBID1
# todo: call in pipeline, specifying which step to start

