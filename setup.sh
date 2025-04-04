
wget -O src/classify_repeat_samps.R https://bitbucket.org/jcignacio/wheat-lab-get-marker-data/raw/d2439277db9257bbb1c8516e637bf4109caa85bc/src/classify_repeat_samps.R
wget -O src/remove_repeats_from_vcf.R https://bitbucket.org/jcignacio/wheat-lab-get-marker-data/raw/d2439277db9257bbb1c8516e637bf4109caa85bc/src/remove_repeats_from_vcf.R

pip3 install matplotlib

R
options(ask = FALSE)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SNPRelate")
BiocManager::install("clstutils")
install.packages("stringr")
q()
