#!/bin/bash

#### go through every directory and do stuff
for d in */ ; do
    cd $d
    # save filenames to list
    filenames=("$(ls -p | grep -v /)")

    # determine mode by looking at the reference file names
    for f in $filenames; do
        if [[ "${f%.*}" == *_consref ]]; then
            mode="emp_consref"
        elif [[ "${f%.*}" == *_QSref ]]; then
            mode="emp_QSref"
        elif [[ "${f%.*}" == *_SGSref ]]; then
            mode="sim_SGSref"
        elif [[ "${f%.*}" == *_simref ]]; then
            mode="sim_simref"
        elif [[ "${f%.*}" == *_NGSref ]]; then
            mode="emp_pureNGS"
        fi
    done
    echo "Mode: $mode"

    # setup data folder
    if [ -d "0000_DATA" ]; then
        echo "Data directory already exists, I skip the setup step..."
    else
        # make dir, enter
        mkdir 0000_DATA
        mv ./*\.* 0000_DATA/ ; cd 0000_DATA

        if [[ $AssemblyRef != "NCBI_HIV1_HXB2.fasta" ]];then
            cp ../../NCBI_HIV1_HXB2.fasta .
        fi

        # make config file with static names, copy files
        configvars=("AssemblyRef" "ShiverConfig" "RefAlignment" "Adapters" "Primers" "Prefix" "vngs_prefix" "VPipeConfig" "SangerCorrection" "start_coord" "end_coord")
        echo "# User supplied variables" >> pipeline.conf
        for c in "${configvars[@]}"; do
            grep $c ../../gabm_batch.conf >> pipeline.conf
            if [[ $c != "Prefix" && $c != "vngs_prefix" && $c != "SangerCorrection" && $c != "start_coord" && $c != "end_coord" ]]; then
                fn=($(grep $c ../../gabm_batch.conf | grep -o -P '(?<=").*(?=")'))
                cp ../../$fn .
            fi
        done

        # handle fastq files
        if [[ $mode == emp_* ]]
        then
            fastqgz_filenames=($(ls ./*.fastq.gz 2>/dev/null)) # make an array out of fastq.gz files in the 0000_DATA directory
            if [ ${#fastqgz_filenames[@]} -eq 0 ] # if you have not found compressed input files
            then
                fastq_filenames=($(ls ./*.fastq 2>/dev/null)) # look for uncompressed input files
                if [ ${#fastq_filenames[@]} -gt 0 ] # if you have found input files
                then 
                    echo "I found the .fastq files"
                    for i in "${fastq_filenames[@]}" # go through every file
                    do
                        gzip $i # compress them
                    done
                    echo "Your .fastq files are now compressed!"
                else
                    echo "There were no .fastq files in the directory!"
                    echo "Exiting..."
                    exit 0
                fi
            else
                echo "I found the .fastq.gz files"
            fi
        fi

        # automatized reference file name setup in pipeline.conf
        echo "# Automatized reference filename variables" >> pipeline.conf
        if [[ $mode == "sim_simref" ]]; then
            ref=("$(find . -name *_simref.fasta)")
            ref=${ref:2}
            echo Master=\"$ref\" >> pipeline.conf
            # simulation parameters
            configvars=("SantaConf" "SimCov")
            echo "# Simulation variables" >> pipeline.conf
            for c in "${configvars[@]}"; do
                echo $c
                grep $c ../../gabm_batch.conf >> pipeline.conf
                if [[ $c != "SimCov" ]]; then
                    fn=($(grep $c ../../gabm_batch.conf | grep -o -P '(?<=").*(?=")'))
                    cp ../../$fn .
                fi
            done
        elif [[ $mode == "sim_SGSref" ]]; then
            ref=("$(find . -name *_SGSref.fasta)")
            ref=${ref:2}
            echo SeqSet=\"$ref\" >> pipeline.conf   
            # simulation parameters
            configvars=("SimCov")
            echo "# Simulation variables" >> pipeline.conf
            for c in "${configvars[@]}"; do
                grep $c ../../gabm_batch.conf >> pipeline.conf
            done     
        elif [[ $mode == "emp_consref" ]]; then
            ref=("$(find . -name *_consref.fasta)")
            ref=${ref:2}
            echo GoldRef=\"$ref\" >> pipeline.conf
        elif [[ $mode == "emp_QSref" ]]; then
            ref=("$(find . -name *_QSref.fasta)")
            ref=${ref:2}
            echo SeqSet=\"$ref\" >> pipeline.conf
        elif [[ $mode == "emp_pureNGS" ]]; then
            ref=("$(find . -name *_NGSref.fasta)")
            ref=${ref:2}
            echo GoldRef=\"$ref\" >> pipeline.conf
        fi

        # automatized read file name setup in pipeline.conf
        echo "# Automatized read filename variables" >> pipeline.conf
        # raw empirical reads
        if [[ $mode == emp_* ]]; then
            fr=("$(find . -name *_R1.fastq.gz)")
            fr=${fr:2}
            echo RawForward=\"$fr\" >> pipeline.conf
            rr=("$(find . -name *_R2.fastq.gz)")
            rr=${rr:2}
            echo RawReverse=\"$rr\" >> pipeline.conf
        fi

        # exit 0000_DATA
        cd ..
    fi

    ################################# run analysis
    bash /home/hazaihiv/Analyses/myscripts/gabm.sh -m $mode |& tee BM_log.txt

    # exit sample dir
    cd ..
done

