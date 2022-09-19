#!/bin/bash
for VARIABLE in L1 L2 L3 S1 S1L1 S1L2 # S2 S2L3
do
    echo ${VARIABLE}
    tail -25 logs/caller_comparison_vaquita_lr/${VARIABLE}_output.log | grep thread
    tail -25 logs/caller_comparison_vaquita_lr/${VARIABLE}_output.log | grep CPU
    tail -25 logs/caller_comparison_vaquita_lr/${VARIABLE}_output.log | grep Elapsed
    tail -25 logs/caller_comparison_vaquita_lr/${VARIABLE}_output.log | grep "Maximum resident"
done
