#!/bin/bash

# Set source and destination directories
src_dir="/home/hazaihiv/Analyses/Benchmarking/RESULTS/batch_runs/SIM/"
dest_dir="/home/hazaihiv/Analyses/Benchmarking/github_repo/data_and_config/SIM_simulated_read/"

# Loop over each directory in the source
for d in "$src_dir"*/ ; do
    # Get the directory name without the full path
    dir_name=$(basename "$d")

    # Create the corresponding directory in the destination
    mkdir -p "$dest_dir/$dir_name"

    # Copy the unique consref_fixed.fasta file
    cp "$d/replicate1/0000_DATA/"*_QSsim_R1.fastq.gz "$dest_dir/$dir_name/" 2>/dev/null
    cp "$d/replicate1/0000_DATA/"*_QSsim_R2.fastq.gz "$dest_dir/$dir_name/" 2>/dev/null

done

echo "File copy completed."

exit 0
