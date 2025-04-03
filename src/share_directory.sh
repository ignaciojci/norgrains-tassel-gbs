#!/bin/bash

source /fs/ess/PAS2444/norgrains-tassel-gbs/config.sh
cd "$WD"

mkdir -p share

ln -s "$WD/raw_VCF" "$WD/share/raw_vcf"
ln -s "$WD/filt" "$WD/share/filtered"
ln -s "$WD/keyfile" "$WD/share/keyfile"

"$RCLONE" copy --copy-links "$WD/share" osu-wheat-lab:norgrains/marker-data/2025_IL
