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
        Rscript ~/Analyses/myscripts/gather_results.R $S $R
        cd ../..
    done
    cd ..
done

# Merge files into 1
Rscript ~/Analyses/myscripts/merge_results.R 
# remove temprary files
shopt -s extglob
rm -v !("GABM_results.csv")