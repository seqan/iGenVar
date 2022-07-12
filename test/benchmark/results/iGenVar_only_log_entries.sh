#!/bin/bash
for VARIABLE in L1 L2 L3 S1 S2 S1L1 S1L2 S1L3 S2L1 S2L2 S2L3
do
    echo ${VARIABLE}
    tail -25 logs/caller_comparison_iGenVar_only/${VARIABLE}_output.log | grep threads
    tail -25 logs/caller_comparison_iGenVar_only/${VARIABLE}_output.log | grep CPU
    tail -25 logs/caller_comparison_iGenVar_only/${VARIABLE}_output.log | grep Elapsed
    tail -25 logs/caller_comparison_iGenVar_only/${VARIABLE}_output.log | grep "Maximum resident"
done
