#!/bin/bash

source /fs/ess/PAS2444/norgrains-tassel-gbs/config.sh
cd "$WD"

mkdir -p share
cd share
mkdir -p raw_VCF filtered keyfile
cd ..

ln -s $WD/raw_VCF/25_norgrains_production.vcf.gz share/raw_vcf/
ln -s $WD/raw_VCF/25_norgrains_production.vcf.gz.csi share/raw_vcf/
ln -s $WD/raw_VCF/25NOR01-09_LIBPREPID_03-26-2025_ReadsPerSample.log share/raw_vcf/

ln -s $WD/filt/05_all_regions_filt.vcf.gz share/filtered/25NOR_all.vcf.gz
ln -s $WD/filt/05_all_regions_filt.vcf.gz.csi share/filtered/25NOR_all.vcf.gz.csi

ln -s $WD/filt/25NORGRAINS_IN.vcf.gz share/filtered_vcf/25NOR_IN.vcf.gz
ln -s $WD/filt/25NORGRAINS_IL.vcf.gz share/filtered_vcf/25NOR_IL.vcf.gz
ln -s $WD/filt/25NORGRAINS_KY.vcf.gz share/filtered_vcf/25NOR_KY.vcf.gz
ln -s $WD/filt/25NORGRAINS_OH.vcf.gz share/filtered_vcf/25NOR_OH.vcf.gz
ln -s $WD/filt/25NORGRAINS_NY.vcf.gz share/filtered_vcf/25NOR_NY.vcf.gz

ln -s $WD/filt/25NORGRAINS_IN.vcf.gz.csi share/filtered_vcf/25NOR_IN.vcf.gz.csi
ln -s $WD/filt/25NORGRAINS_IL.vcf.gz.csi share/filtered_vcf/25NOR_IL.vcf.gz.csi
ln -s $WD/filt/25NORGRAINS_KY.vcf.gz.csi share/filtered_vcf/25NOR_KY.vcf.gz.csi
ln -s $WD/filt/25NORGRAINS_OH.vcf.gz.csi share/filtered_vcf/25NOR_OH.vcf.gz.csi
ln -s $WD/filt/25NORGRAINS_NY.vcf.gz.csi share/filtered_vcf/25NOR_NY.vcf.gz.csi

ln -s $WD/filt/all_samples.txt share/filtered/
ln -s $WD/keyfile/25NOR01-09_LIBPREPID_03-26-2025.txt share/keyfile/

~/softwares/rclone-v1.69.1-linux-amd64/rclone copy --copy-links ./ osu-wheat-lab:norgrains/marker-data/2025/