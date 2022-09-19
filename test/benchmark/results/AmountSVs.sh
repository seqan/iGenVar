#!/bin/bash

echo iGenVar

for input_combination in L1 L2 L3 S1 S2 S1L1 S1L2 S1L3 S2L1 S2L2 S2L3
do
    echo "---------------------------------------"
    case $input_combination in
    S1L1 | S2L1)
        MINQUAL=9
        ;;
    L1)
        MINQUAL=8
        ;;
    S1L2)
        MINQUAL=6
        ;;
    L2 | S2L2)
        MINQUAL=5
        ;;
    S1L3 | S2L3)
        MINQUAL=4
        ;;
    L3)
        MINQUAL=3
        ;;
    S1 | S2)
        MINQUAL=2
        ;;
    *)
        echo "Wrong input combination: $input_combination!"
        ;;
    esac

    echo "${input_combination} - ${MINQUAL}"
    echo "# SV ="
    less results/caller_comparison_iGenVar_only/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "SVTYPE=" | wc -l
    echo "# INS ="
    less results/caller_comparison_iGenVar_only/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "<INS>" | wc -l
    echo "# DUP:TANDEM ="
    less results/caller_comparison_iGenVar_only/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "<DUP:TANDEM>" | wc -l
    echo "# DEL ="
    less results/caller_comparison_iGenVar_only/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "<DEL>" | wc -l
    echo "# INV ="
    less results/caller_comparison_iGenVar_only/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "<INV>" | wc -l
done

echo "------------------------------------------------------------------------------"
echo Vaquita-LR

for input_combination in L1 L2 L3 S1 S1L1 S1L2
do
    echo "---------------------------------------"
    case $input_combination in
    S1 | S1L1 | S1L2)
        MINQUAL=99
        ;;
    L1)
        MINQUAL=97
        ;;
    L2)
        MINQUAL=65
        ;;
    L3)
        MINQUAL=50
        ;;
    *)
        echo "Wrong input combination: $input_combination!"
        ;;
    esac

    echo "${input_combination} - ${MINQUAL}"
    echo "# SV ="
    less results/caller_comparison_vaquita_lr/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "SVTYPE=" | wc -l
    echo "# INS ="
    less results/caller_comparison_vaquita_lr/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "<INS>" | wc -l
    echo "# DUP:TANDEM ="
    less results/caller_comparison_vaquita_lr/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "<DUP:TANDEM>" | wc -l
    echo "# DUP:INT ="
    less results/caller_comparison_vaquita_lr/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "<DUP>" | wc -l
    echo "# DEL ="
    less results/caller_comparison_vaquita_lr/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "<DEL>" | wc -l
    echo "# INV ="
    less results/caller_comparison_vaquita_lr/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "<INV>" | wc -l
    echo "# TRA ="
    less results/caller_comparison_vaquita_lr/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "<TRA>" | wc -l
    echo "# BND ="
    less results/caller_comparison_vaquita_lr/${input_combination}/variants.min_qual_${MINQUAL}.vcf | grep "SVTYPE=BND" | wc -l
done

# ### SVIM:
# `ls -lah results/caller_comparison_long_read/MtSinai_PacBio/SVIM/variants.vcf && ls -lah results/caller_comparison_long_read/MtSinai_PacBio/SVIM/variants.min_qual_1.vcf && ls -lah results/caller_comparison_long_read/PacBio_CCS/SVIM/variants.vcf && ls -lah results/caller_comparison_long_read/PacBio_CCS/SVIM/variants.min_qual_1.vcf && ls -lah results/caller_comparison_long_read/10X_Genomics/SVIM/variants.vcf && ls -lah results/caller_comparison_long_read/10X_Genomics/SVIM/variants.min_qual_1.vcf`
# `grep -o "svim.INS" results/caller_comparison_long_read/.../SVIM/variants.vcf | wc -l`
# `grep -o "svim.DEL" results/caller_comparison_long_read/.../SVIM/variants.vcf | wc -l`
# `grep -o "svim.DUP" results/caller_comparison_long_read/.../SVIM/variants.vcf | wc -l`
# `grep -o "svim.DUP:TANDEM" results/caller_comparison_long_read/.../SVIM/variants.vcf | wc -l`
# `grep -o "svim.DUP:INT" results/caller_comparison_long_read/.../SVIM/variants.vcf | wc -l`
# `grep -o "svim.INV" results/caller_comparison_long_read/.../SVIM/variants.vcf | wc -l`
# `grep -o "svim.BND" results/caller_comparison_long_read/.../SVIM/variants.vcf | wc -l`
# ### Sniffles:
# `ls -lah results/caller_comparison_long_read/MtSinai_PacBio/Sniffles/raw_variants_1.*vcf && ls -lah results/caller_comparison_long_read/PacBio_CCS/Sniffles/raw_variants_1.*vcf && ls -lah results/caller_comparison_long_read/10X_Genomics/Sniffles/raw_variants_1.*vcf`
# `grep -o "SVTYPE=INS" results/caller_comparison_long_read/.../Sniffles/variants.vcf | wc -l`
# `grep -o "SVTYPE=DEL" results/caller_comparison_long_read/.../Sniffles/variants.vcf | wc -l`
# `grep -o "SVTYPE=DUP" results/caller_comparison_long_read/.../Sniffles/variants.vcf | wc -l`
# `grep -o "SVTYPE=INV" results/caller_comparison_long_read/.../Sniffles/variants.vcf | wc -l`
# `grep -o "SVTYPE=INVDUP" results/caller_comparison_long_read/.../Sniffles/variants.vcf | wc -l`
# `grep -o "SVTYPE=TRA" results/caller_comparison_long_read/.../Sniffles/variants.vcf | wc -l`
# ### pbsv & pbsv_without_DUP:
# `grep -o "pbsv.INS" results/caller_comparison_long_read/.../pbsv.../variants.min_qual_1.vcf | wc -l`
# `grep -o "pbsv.DUP" results/caller_comparison_long_read/.../pbsv.../variants.min_qual_1.vcf | wc -l`
# `grep -o "pbsv.INS.DUP" results/caller_comparison_long_read/.../pbsv.../variants.min_qual_1.vcf | wc -l`
# `grep -o "pbsv.INV" results/caller_comparison_long_read/.../pbsv.../variants.min_qual_1.vcf | wc -l`
# `grep -o "pbsv.CNV" results/caller_comparison_long_read/.../pbsv.../variants.min_qual_1.vcf | wc -l`
# ### references:
# `less data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | grep ">" | wc -l`
# 195
# `less data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa | grep ">" | wc -l`
# 3366
# `less data/reference/hg19.fa | grep ">" | wc -l`
# 93
# `less data/reference/hg38.fa | grep ">" | wc -l`
# 455
# `less data/reference/hs37d5.fa | grep ">" | wc -l`
# 86
