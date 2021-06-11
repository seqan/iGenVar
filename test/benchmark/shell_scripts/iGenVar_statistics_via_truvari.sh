#!/usr/bin/env sh
current_date=$(date +"%Y-%m-%d")
echo "Logifle written to: ${current_date}_statistics_via_truvari.log"


# -------- -------- Prepare logfile and runtime computation -------- -------- #

exec 3>&1 4>&2 # store original streams
trap 'exec 2>&4 1>&3' 0 1 2 3 # restore original streams when scripts ends
exec 1>logs/${current_date}_statistics_via_truvari.log 2>&1 # write stdout and stdcerr to logfile

set -ex

echo $(date)
start=$(date +%s) # get starting date

# -------- -------- pre installation steps -------- -------- #
# for truvari you need python3,
# we recommend to use miniconda:
# ------------------------------
# wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
# chmod +x Miniconda2-latest-Linux-x86_64.sh
# sh Miniconda2-latest-Linux-x86_64.sh

# conda env create -f environment.yml

# source ~/.bashrc                                    # if the .bashrc is not executed automatically
conda activate iGenVar_benchmark

# -------- -------- get truth set ../../data -------- -------- #
cd ../../data
mkdir -p truth_set && cd truth_set

wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/../../data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/../../data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose \
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/../../data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed

echo "$(tput setaf 1)$(tput setab 7)------- truth set downloaded (1/3) --------$(tput sgr 0)" 1>&3

cd ../..


# -------- -------- init for truvari -------- -------- #

cd results/${current_date}

for min_var_length in 10 30 50 100 500 1000 2000
do
    echo "Looping min_var_length number $min_var_length"

    vcf_file_name=min_var_length_${min_var_length}_output.vcf
    bgzip -c ${vcf_file_name} > ${vcf_file_name}.gz
    tabix -p vcf ${vcf_file_name}.gz
done


for max_var_length in 500000 1000000 2000000
do
    echo "Looping max_var_length number $max_var_length"

    vcf_file_name=max_var_length${max_var_length}_output.vcf
    bgzip -c ${vcf_file_name} > ${vcf_file_name}.gz
    tabix -p vcf ${vcf_file_name}.gz
done

for max_tol_inserted_length in 1 5 10
do
    echo "Looping max_tol_inserted_length number $max_tol_inserted_length"

    vcf_file_name=max_tol_inserted_length_${max_tol_inserted_length}_output.vcf
    bgzip -c ${vcf_file_name} > ${vcf_file_name}.gz
    tabix -p vcf ${vcf_file_name}.gz
done

for max_overlap in 5 10 20
do
    echo "Looping max_overlap number $max_overlap"

    vcf_file_name=max_overlap_${max_overlap}_output.vcf
    bgzip -c ${vcf_file_name} > ${vcf_file_name}.gz
    tabix -p vcf ${vcf_file_name}.gz
done

for hierarchical_clustering_cutoff in 10 50 100 150
do
    echo "Looping hierarchical_clustering_cutoff number $hierarchical_clustering_cutoff"

    vcf_file_name=hc_cutoff_${hierarchical_clustering_cutoff}_output.vcf
    bgzip -c ${vcf_file_name} > ${vcf_file_name}.gz
    tabix -p vcf ${vcf_file_name}.gz
done

echo "$(tput setaf 1)$(tput setab 7)------- truvari init done (2/3) --------$(tput sgr 0)" 1>&3

# -------- -------- run truvari -------- -------- #

for min_var_length in 10 30 50 100 500 1000 2000
do
    echo "Looping min_var_length number $min_var_length"
    vcf_file_name=min_var_length_${min_var_length}_output.vcf

    truvari bench -b ../../data/truth_set/HG002_SVs_Tier1_v0.6.vcf.gz -c ${vcf_file_name}.gz \
        -o min_var_length_${min_var_length}_truvari_default \
        --includebed ../../data/truth_set/HG002_SVs_Tier1_v0.6.bed -p 0

    truvari bench -b ../../data/truth_set/HG002_SVs_Tier1_v0.6.vcf.gz -c ${vcf_file_name}.gz \
        -o min_var_length_${min_var_length}_truvari_multimatch \
        --includebed ../../data/truth_set/HG002_SVs_Tier1_v0.6.bed -p 0 --multimatch
done


for max_var_length in 500000 1000000 2000000
do
    echo "Looping max_var_length number $max_var_length"
    vcf_file_name=max_var_length_${max_var_length}_output.vcf

    truvari bench -b ../../data/truth_set/HG002_SVs_Tier1_v0.6.vcf.gz -c ${vcf_file_name}.gz \
        -o max_var_length_${max_var_length}_truvari_default \
        --includebed ../../data/truth_set/HG002_SVs_Tier1_v0.6.bed -p 0

    truvari bench -b ../../data/truth_set/HG002_SVs_Tier1_v0.6.vcf.gz -c ${vcf_file_name}.gz \
        -o max_var_length_${max_var_length}_truvari_multimatch \
        --includebed ../../data/truth_set/HG002_SVs_Tier1_v0.6.bed -p 0 --multimatch
done

for max_tol_inserted_length in 1 5 10
do
    echo "Looping max_tol_inserted_length number $max_tol_inserted_length"
    vcf_file_name=max_tol_inserted_length_${max_tol_inserted_length}_output.vcf

    truvari bench -b ../../data/truth_set/HG002_SVs_Tier1_v0.6.vcf.gz -c ${vcf_file_name}.gz \
        -o max_tol_inserted_length_${max_tol_inserted_length}_truvari_default \
        --includebed ../../data/truth_set/HG002_SVs_Tier1_v0.6.bed -p 0

    truvari bench -b ../../data/truth_set/HG002_SVs_Tier1_v0.6.vcf.gz -c ${vcf_file_name}.gz \
        -o max_tol_inserted_length_${max_tol_inserted_length}_truvari_multimatch \
        --includebed ../../data/truth_set/HG002_SVs_Tier1_v0.6.bed -p 0 --multimatch
done

for max_overlap in 5 10 20
do
    echo "Looping max_overlap number $max_overlap"
    vcf_file_name=max_overlap_${max_overlap}_output.vcf

    truvari bench -b ../../data/truth_set/HG002_SVs_Tier1_v0.6.vcf.gz -c ${vcf_file_name}.gz \
        -o max_overlap_${max_overlap}_truvari_default \
        --includebed ../../data/truth_set/HG002_SVs_Tier1_v0.6.bed -p 0

    truvari bench -b ../../data/truth_set/HG002_SVs_Tier1_v0.6.vcf.gz -c ${vcf_file_name}.gz \
        -o max_overlap_${max_overlap}_truvari_multimatch \
        --includebed ../../data/truth_set/HG002_SVs_Tier1_v0.6.bed -p 0 --multimatch
done

for hierarchical_clustering_cutoff in 10 50 100 150
do
    echo "Looping hierarchical_clustering_cutoff number $hierarchical_clustering_cutoff"
    vcf_file_name=hc_cutoff_${hierarchical_clustering_cutoff}_output.vcf

    # find ./ -type d -name "hc_cutoff_${hierarchical_clustering_cutoff}_truvari_default" -exec rm -rf {} +        # delete if already exists
    truvari bench -b ../../data/truth_set/HG002_SVs_Tier1_v0.6.vcf.gz -c ${vcf_file_name}.gz \
        -o hc_cutoff_${hierarchical_clustering_cutoff}_truvari_default \
        --includebed ../../data/truth_set/HG002_SVs_Tier1_v0.6.bed -p 0

    # find ./ -type d -name "hc_cutoff_${hierarchical_clustering_cutoff}_truvari_multimatch" -exec rm -rf {} +     # delete if already exists
    truvari bench -b ../../data/truth_set/HG002_SVs_Tier1_v0.6.vcf.gz -c ${vcf_file_name}.gz \
        -o hc_cutoff_${hierarchical_clustering_cutoff}_truvari_multimatch \
        --includebed ../../data/truth_set/HG002_SVs_Tier1_v0.6.bed -p 0 --multimatch
done

cd ../..

echo "$(tput setaf 1)$(tput setab 7)------- truvari default and multimatch statistics calculated (3/3) --------$(tput sgr 0)" 1>&3


# ----------------------------------------
echo "$(tput setaf 1)$(tput setab 7)------- Statistics calculated via truvari - done --------$(tput sgr 0)" 1>&3

end=$(date +%s)                 # get end-date
runtime=$(((end-start)/60))     # calculate runtime
echo "$runtime minutes"
echo "$(tput setaf 1)$(tput setab 7)This run took $runtime minutes$(tput sgr 0)" 1>&3
