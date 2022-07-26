#!/bin/bash
for VARIABLE in L1 L2 L3 S1 S2 S1L1 S1L2 S1L3 S2L1 S2L2 S2L3
do
    echo "---------------------------------------"
    echo ${VARIABLE}
    for MINQUAL in 1 2 3 4 5 6 7 8 9
    do
        echo ${MINQUAL}
        less results/caller_comparison_iGenVar_only/eval/${VARIABLE}/no_DUP_and_INV.min_qual_${MINQUAL}/summary.txt | grep f1
        less results/caller_comparison_iGenVar_only/eval/${VARIABLE}/no_DUP_and_INV.min_qual_${MINQUAL}/summary.txt | grep precision
        less results/caller_comparison_iGenVar_only/eval/${VARIABLE}/no_DUP_and_INV.min_qual_${MINQUAL}/summary.txt | grep recall
        less results/caller_comparison_iGenVar_only/${VARIABLE}/variants.no_DUP_and_INV.min_qual_${MINQUAL}.vcf | grep "SVTYPE=" | wc -l
    done
done
