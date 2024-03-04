#!/bin/bash

#### go through every directory and generate result table for the scenario
for S in */ ; do
    cd $S
    for R in */ ; do
        cd $R
        cd 1000_BMresults 
        #modify _tm.txt files
        sed -i 's/:/;/2g' *_tm.txt
        #gather resuts
        Rscript ~/Analyses/myscripts/gather_results.R $S $R &>> ../../../gather_results.txt
        cd ../..
    done
    cd ..
done

# Merge files into 1
Rscript ~/Analyses/myscripts/merge_results.R &> merge_results.txt
# remove temprary files
shopt -s extglob
rm -v !("GABM_NGS_results.csv"|"GABM_SSNGS_results.csv"|"GABM_results.csv"|"GABM_SGS_results.csv"|"gather_results.txt"|"merge_results.txt")