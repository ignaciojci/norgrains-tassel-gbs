#!/bin/bash

## Filter VCF file - parallel
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script filters a VCF file variant-wise. It will also subset samples based
## on a user-defined list of sample names, if supplied.
##
## The script takes two positional, command-line arguments:
##   1) The path to a .bed file defining regions of VCF file to extract and filter.
##      This file must be sorted by chromosome and then start position, and should
##      not contain any overlapping regions.
##   2) An integer, defining which line of the .bed file to use. Supplying 0 will
##      use all lines of the .bed file. This integer will typically be supplied
##      by a parallel dispatch script, which will allow for filtering multiple
##      regions simultaneously.
##
## You can create a .bed file of whole chromosomes/contigs from the contig lines in the VCF/BCF:
##
##   bcftools view -h <file.bcf> | grep "contig=" | sed -e 's/^.*ID=//' -e 's/,length=/\t0\t/' -e 's/>//' > file.bed
##
## All other constants are set within the script, below this description.
##
## Whether the input is a VCF or BCF file, the output will always be a BCF file.
## If using all regions, the output file will be placed in the user-defined
## output directory as "all_regions.bcf". If a single region is filtered, then
## the output BCF file will have a name formatted as [region_integer]_[chrom]_[start]_[end].bcf
## The integer in the file name is formatted with leading zeroes to allow concatenation
## of multiple BCF files without a sorting step. For instance, an output BCF file name
## of region number 15 out of 1000 might be formatted as:
##
##   0015_2A_100000_200000.bcf
##
## NOTES:
##
##   The max_miss parameter works in an opposite manner of the VCFTools
##   implementation, which I find unintuitive. Here, if max_miss is set to 1, SNPs with 100%
##   missing data would theoretically be allowed, while setting max_miss to 0 will only
##   allow SNPs without any missing data
##
##   Samples in the output file will always be sorted.
##
##   Filtering to only retain biallelic variants, only retain SNPs, etc. must
##   be hand-coded into the first bcftools call in the piped command
##
##   The input VCF/BCF file must be indexed, and must include contig lines in
##   the header.
##
##   A regions .bed file is still required if filtering an entire VCF. This can
##   be created from the contig lines in the VCF header, e.g.:
##
##     bcftools view -h [input.vcf.gz] | grep "##contig" | sed -e 's/##contig=<ID=//' -e 's/,length=/\t0\t/' -e 's/[,>$]//'
##
##   For depth filtering, this script uses the DP value in the INFO column. This
##   value is the total depth of reads that passed the quality control parameters
##   defined during SNP calling, (e.g. minimum mapq value). FORMAT DP values, if
##   present, are for raw read counts per sample.
##
##   snpgap will remove SNPs within n bases of an indel (or overlapping)
##   indelgap will thin clusters of indels within n bases of each other, to only retain one
##
##   BCFTools is designed to work with .bcf files. If .vcf files are supplied, it
##   first converts them to .bcf format internally. Therefore, using a .bcf file
##   as input is faster, though BCF files may be larger than .vcf.gz files.
##
##   The output file will contain some additional INFO tags, including MAF, F_MISSING
##   (the frequency of missing data), AC_Het (The count of reference alleles in
##    heterozygous calls), and NS (number of samples with data)
##
##   The presence of multi-allelic SNPs may make the het. filtering slightly inaccurate.
##   The current implementation relies on a het. call always containing one REF
##   allele. However, if there is more than one ALT allele, then a variant may
##   contain 0 REF alleles while still being het.
################################################################################


#### SLURM job control ####

#SBATCH --time=10:00:00
#SBATCH --nodes=1 --ntasks-per-node=20 --partition=hugemem
#SBATCH --job-name=TestProduction0404
#SBATCH --account=PAS2444
#SBATCH --mail-user=hannahmaohuang@gmail.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error

set -u -e -o pipefail -x

#### User-defined constants ####################################################

## Note that SNP depth and proportion of missing data are highly correlated

vcf_inp=$1
out_dir=$2
samp_file="none"
min_maf=0.05
max_miss=0.2
max_het=0.25
min_mean_nonmiss_dp=0
max_mean_nonmiss_dp=10000
het2miss="false"
snpgap=3
indelgap=3


#### Executable ################################################################

## Using conda here to utilize bcftools plugins and plot-vcfstats
#module load singularity/3.7.1
#module load bcftools
#module load bcftools
#module load miniconda
#source activate htslib

echo
echo "Start time:"
date

regions_bed=$3
array_ind=$4
nthreads=$5

## Get first letter of true/false het2miss option and
## Set some constants depending on het2miss
het2miss="${het2miss:0:1}"
if [[ "$het2miss" == [Tt] ]]; then
	max_het="NA"
	het_string="./."
fi

## Output filtering parameters and generate the samples file
if [[ $array_ind -eq 0 || $array_ind -eq 1 ]]; then
    mkdir -p "$out_dir"
    mkdir -p "${out_dir}/temp_files"

    ## Echo input parameters to output dir
    echo -e "Input VCF\t${vcf_inp}" > "${out_dir}/varwise_filt_params.txt"
    echo -e "Output directory\t${out_dir}" >> "${out_dir}/varwise_filt_params.txt"
    echo -e "Regions .bed file\t${regions_bed}" >> "${out_dir}/varwise_filt_params.txt"
    echo -e "Sample subset list\t${samp_file}" >> "${out_dir}/varwise_filt_params.txt"
    echo -e "Minimum MAF\t${min_maf}" >> "${out_dir}/varwise_filt_params.txt"
    echo -e "Max missing proportion\t${max_miss}" >> "${out_dir}/varwise_filt_params.txt"
    echo -e "Max het. proportion\t${max_het}" >> "${out_dir}/varwise_filt_params.txt"
    echo -e "Min. mean depth in samples with data\t${min_mean_nonmiss_dp}" >> "${out_dir}/varwise_filt_params.txt"
    echo -e "Max mean depth in samples with data\t${max_mean_nonmiss_dp}" >> "${out_dir}/varwise_filt_params.txt"
    echo -e "Heterozygous calls converted to missing data?\t${het2miss}" >> "${out_dir}/varwise_filt_params.txt"
    echo -e "SNP overlap gap\t${snpgap}" >> "${out_dir}/varwise_filt_params.txt"
    echo -e "Indel overlap gap\t${indelgap}" >> "${out_dir}/varwise_filt_params.txt"

    ## If samp_file exists, use it to subset samples
    ## Otherwise retain all samples present in the input file
    if [[ -f "$samp_file" ]]; then
        cut -f 1 "$samp_file" | sort > "${out_dir}/temp_files/samp_file.txt"
    else
        bcftools query --list-samples "$vcf_inp" | sort > "${out_dir}/temp_files/samp_file.txt"
    fi
else
    sleep 10s
fi

## If the supplied array_ind is 0, we copy full regions .bed file to temp dir.
## Otherwise extract single region from .bed file
if [[ array_ind -eq 0 ]]; then
    label="05_all_regions_filt"
    cp "$regions_bed" "${out_dir}/temp_files/${label}.bed"
else
    ## This ensures output VCFs will have names in order, to avoid needing to
    ## sort after concattenating them together
    n=$(($(wc -l < "$regions_bed" | wc -c) - 1))
    prefix=$(printf "%0${n}d" $array_ind)

    region=$(head -n $array_ind "$regions_bed" | tail -n 1)
    suffix=$(echo "$region" | tr "\t" "_")
    label="${prefix}_${suffix}"
    echo "$region" > "${out_dir}/temp_files/${label}.bed"
fi



## Big copy-pasted if-else block for different actions depending het2miss
## Using het2miss will leave you with some SNPs that only have one allele, so some extra steps
## Are necessary to get rid of these
echo "Filtering VCF..."
echo

if [[ "$het2miss" == [Tt] ]]; then
    bcftools view "${vcf_inp}" \
        --samples-file "${out_dir}/temp_files/samp_file.txt" \
        --min-alleles 2 \
        --max-alleles 2 \
        --types snps \
        --regions-file "${out_dir}/temp_files/${label}.bed" \
        --output-type u |
    bcftools +setGT --output-type u - -- --target-gt q \
        --include 'GT="het"' \
        --new-gt "$het_string" |
    bcftools +fill-tags --output-type u - -- -t MAF,F_MISSING,AC_Het,NS |
    bcftools +fill-tags --output-type u - -- -t 'DPsum:1=int(sum(FORMAT/DP))' |
    bcftools view - \
        --exclude "INFO/F_MISSING > ${max_miss} || INFO/MAF < ${min_maf} || (INFO/DPsum)/(INFO/NS) < ${min_mean_nonmiss_dp} || (INFO/DPsum)/(INFO/NS) > ${max_mean_nonmiss_dp} || (INFO/AC_Het)/(INFO/NS) > ${max_het}" \
        --output-type u |
    bcftools filter - \
        --SnpGap $snpgap \
        --IndelGap $indelgap \
        --output-type z \
        --threads $nthreads \
        --output "${out_dir}/${label}.vcf.gz"
else
    bcftools view "${vcf_inp}" \
        --samples-file "${out_dir}/temp_files/samp_file.txt" \
        --min-alleles 2 \
        --max-alleles 2 \
        --types snps \
        --regions-file "${out_dir}/temp_files/${label}.bed" \
        --output-type u |
    bcftools +fill-tags --output-type u - -- -t MAF,F_MISSING,AC_Het,NS |
    bcftools +fill-tags --output-type u - -- -t 'DPsum:1=int(sum(FORMAT/DP))' |
    bcftools view - \
        --exclude "INFO/F_MISSING > ${max_miss} || INFO/MAF < ${min_maf} || (INFO/DPsum)/(INFO/NS) < ${min_mean_nonmiss_dp} || (INFO/DPsum)/(INFO/NS) > ${max_mean_nonmiss_dp} || (INFO/AC_Het)/(INFO/NS) > ${max_het}" \
        --output-type u |
    bcftools filter - \
        --SnpGap $snpgap \
        --IndelGap $indelgap \
        --output-type z \
        --threads $nthreads \
        --output "${out_dir}/${label}.vcf.gz"
fi

## If a whole VCF/BCF file was filtered, generate some summary stats and plots
if [[ $array_ind -eq 0 ]]; then
    bcftools index -c "${out_dir}/${label}.vcf.gz"

    ## Generate summary stats
    bcftools stats "${out_dir}/${label}.vcf.gz" > "${out_dir}/stats.txt"

    ## Create plots from summary stats
    mkdir "${out_dir}/plots"
    plot-vcfstats --prefix "${out_dir}/plots" --no-PDF "${out_dir}/stats.txt"
fi

## Generate summary stats using TASSEL
#echo "Generating summary statistics with TASSEL..."
#$TASSEL_PL -vcf "${out_dir}/${label}.vcf.gz" \
#           -GenotypeSummaryPlugin \
#           -endPlugin \
#           -export "${out_dir}"/summary

echo
echo "End time:"
date
