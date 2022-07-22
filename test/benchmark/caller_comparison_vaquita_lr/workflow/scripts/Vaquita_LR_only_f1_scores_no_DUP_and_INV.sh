#!/bin/bash
for VARIABLE in L1 L2 L3 S1 S1L1 S1L2 # S2 S2L3
do
    echo "---------------------------------------"
    echo ${VARIABLE}
    for MINQUAL in 49 50 51 64 65 66 96 97 98 99
    do
        echo ${MINQUAL}
        less results/caller_comparison_vaquita_lr/eval/${VARIABLE}/DUP_as_INS.min_qual_${MINQUAL}/summary.txt | grep f1
        less results/caller_comparison_vaquita_lr/eval/${VARIABLE}/DUP_as_INS.min_qual_${MINQUAL}/summary.txt | grep precision
        less results/caller_comparison_vaquita_lr/eval/${VARIABLE}/DUP_as_INS.min_qual_${MINQUAL}/summary.txt | grep recall
        less results/caller_comparison_vaquita_lr/${VARIABLE}/variants.DUP_as_INS.min_qual_${MINQUAL}.vcf | grep "SVTYPE=" | wc -l
    done
done
