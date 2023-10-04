#!/bin/bash

##############################
# Help                       #
##############################
Help()
{
   # Display Help
    echo "## Benchmarking of HIV-1 genome assembly pipelines using paired-end short-read sequencing data ##"
    echo
    echo "Start this script from the analysis directory. Copy all input files into the 0000_DATA subdirectory."
    echo "Specify file names in the config file (pipeline.conf)"
    echo
    echo "Dependencies:"
    echo "   Command line tools available in package managers: mafft, samtools, bamtools, smalt, tabix, fastqc, trimmomatic, seqtk"
    echo "   Docker images: quay.io/biocontainers/ngshmmalign:0.1.1--ha04c180_4, kboltonlab/lofreq:latest, staphb/quast:latest, biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1, biocontainers/art-nextgen-simulation-tools:v20160605dfsg-3-deb_cv1, fazekasda/shiver:22.04, broadinstitute/viral-ngs:latest"
    echo "   Other command line tools: bbmap, picard, santa-sim, anaconda, smaltalign, V-Pipe"
    echo
    echo "Syntax: path/gabm.sh [-m|h]"
    echo
    echo "options:"
    echo "m     Mode ('emp_QSref' - empirical data with known sequences,'emp_consref' - empirical data with unknown sequences,'sim' - simulated data)"
    echo "h     Print this Help."
    echo
}

####################################
###   Process the input options. ###
####################################

## Get the options
while getopts ":hm:" option; do
    case $option in
        h) # display Help
            Help
            exit;;
        m) # Enter a mode
            Mode=$OPTARG
            Accepted=("emp_QSref" "emp_consref" "sim_simref" "sim_SGSref")
            if [[ " ${Accepted[*]} " =~ " ${Mode} " ]]
            then
                echo "Mode: $Mode"
            else
                echo "Error: Invalid mode ($Mode)"
                exit
            fi
            ;;
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done

## Sourcing config file
. ./0000_DATA/pipeline.conf # source config file containing file names

## Handling .fastq files
if [[ $Mode == "emp_QSref" || $Mode == "emp_consref" ]]
then
    fastqgz_filenames=($(ls 0000_DATA/*.fastq.gz 2>/dev/null)) # make an array out of fastq.gz files in the 0000_DATA directory
    if [ ${#fastqgz_filenames[@]} -eq 0 ] # if you have not found compressed input files
    then
        fastq_filenames=($(ls 0000_DATA/*.fastq 2>/dev/null)) # look for uncompressed input files
        if [ ${#fastq_filenames[@]} -gt 0 ] # if you have found input files
        then 
            for i in "${fastq_filenames[@]}" # go through every file
            do
                gzip $i # compress them
            done
            echo "Your .fastq files are now compressed!"
            RawForward=$RawForward\.gz # rename bash variable
            RawReverse=$RawReverse\.gz
        else
            echo "There were no .fastq files in the directory!"
            echo "Put input files into the ./0000_DATA directory according to the manual. Exiting..."
            exit 0
        fi
    fi
fi

## Handling contamination folder
#### determine contamination scenario
if [ -d "0000_contamination" ]
then
    Contamination=true
else
    Contamination=false
fi


############################
## Benchmarking baseline  ##
############################

## QS sequence generation with santa-sim
if [[ $Mode == "sim_simref" ]]
then
    cd 0000_DATA
    # correct consensus sequences by replacing Ns
    python3 ~/Analyses/myscripts/changeNbases.py $Master
    rm $Master; mv genome_mod.fasta $Master
    # define generation in Santaconf
    Generation=$(shuf -i 50-1500 -n 1)
    sed -i "s#<generationCount>1000</generationCount>#<generationCount>'$Generation'</generationCount>#" $SantaConf
    sed -i "s#<atGeneration>1000</atGeneration>#<atGeneration>'$Generation'</atGeneration>#" $SantaConf
    sed -i "s/>'/>/g" $SantaConf
    sed -i "s/'</</g" $SantaConf
    # simulation
    java -jar /home/hazaihiv/Software/santa-sim/dist/santa.jar $SantaConf
    mv alignment_1.fasta "${Master%.*}_QSsim.fasta"
    cd ..
    SeqSet="${Master%.*}_QSsim.fasta" # define sequence set for the 'sim' mode
fi

## Processing QS sequences
mkdir -p 0010_goldenseqs # make dir
if [[ $Mode != "emp_consref" ]]
then
    # multiple alignment
    # mafft --localpair --maxiterate 1000 ./0000_DATA/"${SeqSet}" > ./0010_goldenseqs/"${SeqSet%.*}"_aligned.fasta
    mafft ./0000_DATA/"${SeqSet}" > ./0010_goldenseqs/"${SeqSet%.*}"_aligned.fasta # just for development purposes
    # consensus
    sudo chmod +777 -R ../* # solve permission issue with output files
    # call consensus from the QS sequence set
    docker run --rm -it -v $(pwd)/:/data biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1 /bin/bash -c "em_cons -sequence ./0010_goldenseqs/'${SeqSet%.*}'_aligned.fasta -outseq ./0010_goldenseqs/'${SeqSet%.*}'_aligned_consensus.fasta"
    samtools faidx ./0010_goldenseqs/"${SeqSet%.*}"_aligned_consensus.fasta # index consensus
fi

## Define golden reference sequence
if [[ $Mode != "emp_consref" ]]
then
    GoldRef="${SeqSet%.*}"_aligned_consensus.fasta # consensus of known sequences
else
    sudo chmod +777 -R ../*
    cp ./0000_DATA/$GoldRef ./0010_goldenseqs # arbitrary sequence (results might indicate outlier assemblies between genome assembly software)
fi
python3 ~/Analyses/myscripts/cropGenomeByMotifs.py ./0010_goldenseqs/$GoldRef ./0000_DATA/$RegionBoundaries ./0010_goldenseqs/"${GoldRef%.*}"_fixed.fasta
GoldRef="${GoldRef%.*}"_fixed.fasta

## Simulation of PE NGS reads
if [[ $Mode == "sim_simref" || $Mode == "sim_SGSref" ]]
then
    # calculate coverage
    SeqNum=("$(grep ">" ./0010_goldenseqs/"${SeqSet%.*}"_aligned.fasta | wc -l)")
    SimPerseqCov=("$(expr $SimCov / $SeqNum)")
    echo "Number of reads per sequence:" $SimPerseqCov
    # generate sequences with ART
	docker run --rm -it -v $(pwd)/:/data 19da3c67bd0edccc378d135e79b7ce8a376d3aabccfb087c97cb340920ff1eed /bin/bash -c "cd 0000_DATA; art_illumina -sam -i '$SeqSet' -l 250 -s 1 -f '$SimPerseqCov' -ss MSv3 -o '${SeqSet%.*}' -m 700 -nf 0 -p -qL 20 -qU 40 -rs 1680244882"
    # remove junk files
    sudo rm ./0000_DATA/"${SeqSet%.*}"1.aln ./0000_DATA/"${SeqSet%.*}"2.aln ./0000_DATA/"${SeqSet%.*}".sam
    # move results to the data directory
    mv ./0000_DATA/"${SeqSet%.*}"1.fq ./0000_DATA/"${SeqSet%.*}"_R1.fastq; mv ./0000_DATA/"${SeqSet%.*}"2.fq ./0000_DATA/"${SeqSet%.*}"_R2.fastq
    # compress the simulated reads
    gzip ./0000_DATA/"${SeqSet%.*}"_R1.fastq; gzip ./0000_DATA/"${SeqSet%.*}"_R2.fastq
    # name read input variables
    RawForward="${SeqSet%.*}"_R1.fastq.gz
    RawReverse="${SeqSet%.*}"_R2.fastq.gz
fi

## Preprocessing PE NGS reads 
# make dir
mkdir -p 0020_procreads
# QC
fastqc -o 0020_procreads ./0000_DATA/$RawForward
fastqc -o 0020_procreads ./0000_DATA/$RawReverse
# Read trimming and filtering with trimmomatic
TrimmomaticPE ./0000_DATA/$RawForward ./0000_DATA/$RawReverse ./0020_procreads/"${RawForward%.*.*}"_paired.fastq.gz ./0020_procreads/"${RawForward%.*.*}"_unpaired.fastq.gz ./0020_procreads/"${RawReverse%.*.*}"_paired.fastq.gz ./0020_procreads/"${RawReverse%.*.*}"_unpaired.fastq.gz ILLUMINACLIP:./0000_DATA/$Adapters:2:10:7:1:true MINLEN:50 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20
gzip ./0020_procreads/"${RawForward%.*.*}"_paired.fastq; gzip ./0020_procreads/"${RawReverse%.*.*}"_paired.fastq
# name preprocessed reads
PreprocessedForward="${RawForward%.*.*}"_paired.fastq.gz
PreprocessedReverse="${RawReverse%.*.*}"_paired.fastq.gz

## Golden mapping
# make dir
mkdir -p 0030_goldenmap
# index consensus and map empirical reads with smalt to the consensus
smalt index -k 15 -s 3 ./0010_goldenseqs/"${GoldRef%.*}" ./0010_goldenseqs/$GoldRef
smalt map -x -y 0.7 -j 0 -i 2000 -o ./0030_goldenmap/"${GoldRef%.*}"_goldenmap.sam ./0010_goldenseqs/"${GoldRef%.*}" ./0020_procreads/$PreprocessedForward ./0020_procreads/$PreprocessedReverse 
# sort and index the bam file with samtools and bamtools
samtools sort ./0030_goldenmap/"${GoldRef%.*}"_goldenmap.sam -o ./0030_goldenmap/"${GoldRef%.*}"_goldenmap_sorted.bam
bamtools index -in ./0030_goldenmap/"${GoldRef%.*}"_goldenmap_sorted.bam
# deduplication with picard, indexing with bamtools
sudo picard-tools MarkDuplicates INPUT=./0030_goldenmap/"${GoldRef%.*}"_goldenmap_sorted.bam OUTPUT=./0030_goldenmap/"${GoldRef%.*}"_goldenmap_sorted_dedupreads.bam METRICS_FILE=./0030_goldenmap/"${GoldRef%.*}"_goldenmap_sorted_dedupmetrics.txt
bamtools index -in ./0030_goldenmap/"${GoldRef%.*}"_goldenmap_sorted_dedupreads.bam

## Golden variant calling
# make dir
mkdir -p 0040_goldenvar
# call variants from this alignment qith lofreq, compress qith bgzip, index with tabix
docker run --rm -it -v $(pwd)/:/data kboltonlab/lofreq /bin/bash -c "cd ../data; lofreq call --call-indels -f ./0010_goldenseqs/'$GoldRef' -o ./0040_goldenvar/'${GoldRef%.*}'.vcf ./0030_goldenmap/'${GoldRef%.*}'_goldenmap_sorted_dedupreads.bam"
bgzip -c ./0040_goldenvar/"${GoldRef%.*}".vcf > ./0040_goldenvar/"${GoldRef%.*}".vcf.gz
tabix -fp vcf ./0040_goldenvar/"${GoldRef%.*}".vcf.gz

## Add contamination to the assembly datasets
if [ $Contamination = true ]
then
    ## Simulation of contaminated read sets
    cd 0000_contamination
    Contamfile=("$(ls -p | grep -v /)")
    docker run --rm -it -v $(pwd)/:/data 19da3c67bd0edccc378d135e79b7ce8a376d3aabccfb087c97cb340920ff1eed /bin/bash -c "art_illumina -sam -i '$Contamfile'  -l 250 -s 1 -c 8600 -ss MSv3 -o '${Contamfile%.*}' -m 700 -nf 0 -p -qL 20 -qU 40 -rs 1680244882"
    gzip "${Contamfile%.*}"1.fq
    gzip "${Contamfile%.*}"2.fq
    ForwardContam="${Contamfile%.*}"1.fq.gz
    ReverseContam="${Contamfile%.*}"2.fq.gz
    cd ..
    # add contaminated reads to the dataset
    cat 0000_DATA/$RawForward 0000_contamination/$ForwardContam > 0000_DATA/"${RawForward%_*.*.*}"_cont_R1.fastq.gz
    RawForwardMod="${RawForward%_*.*.*}"_cont_R1.fastq.gz
    cat 0000_DATA/$RawReverse 0000_contamination/$ReverseContam > 0000_DATA/"${RawReverse%_*.*.*}"_cont_R2.fastq.gz
    RawReverseMod="${RawReverse%_*.*.*}"_cont_R2.fastq.gz
else   
    RawForwardMod=$RawForward
    RawReverseMod=$RawReverse
fi
# interleave forward and reverse reads (for smaltalign):
reformat.sh in1=./0000_DATA/$RawForwardMod in2=./0000_DATA/$RawReverseMod out=./0020_procreads/"${RawForwardMod%_*.*.*}"_R1_R2.fastq.gz


#################################
###      Genome assembly      ###
#################################

# make directories
mkdir -p 0100_shiver
mkdir -p 0200_smaltalign
mkdir -p 0300_viralngs
mkdir -p 0400_vpipe/config

## Shiver
# specify read names to shiver
echo ForwardReads=\"$RawForwardMod\" >> ./0000_DATA/pipeline.conf
echo ReverseReads=\"$RawReverseMod\" >> ./0000_DATA/pipeline.conf
# copy input files to shiver directory
cp ./0000_DATA/$RefAlignment ./0100_shiver # alignment of reference sequences
cp ./0000_DATA/$Adapters ./0100_shiver; cp ./0000_DATA/$Primers ./0100_shiver # adapters and primers
cp ./0000_DATA/$ShiverConfig ./0100_shiver; cp ./0000_DATA/pipeline.conf ./0100_shiver # shiver config files
cp ./0000_DATA/$RawForwardMod ./0100_shiver; cp ./0000_DATA/$RawReverseMod ./0100_shiver # processed input reads (trimmomatic inside the shiver pipeline, but same specification)
# go inside directory
cd 0100_shiver
# run shiver pipeline
docker run --rm -it -v $(pwd)/:/data hcovlab/dshiver:latest full
# remove input files and exit directory
rm $RefAlignment $Adapters $ShiverConfig pipeline.conf $Primers $RawForwardMod $RawReverseMod
cd ..
# create dummy fasta if pipeline fails
if [ ! -f 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30.fasta ]
then
    echo "Shiver genome doesn't exist. Prepare dummy genome..."
    echo ">shiver_error" > 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30_fixed.fasta
    printf 'N%.0s' {1..8500} >> 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30_fixed.fasta
    echo >> 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30_fixed.fasta
else
    echo "Shiver genome exists. Trimming..."
    # shiver consensus file contains 2 sequences, remove one without reference imputation, cut out specified region
    cat 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30.fasta | awk '/>'"${Prefix}"'_remap_consensus/ {getline; while(!/>/) {getline}} 1' > 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30_prefixed.fasta
    python3 ~/Analyses/myscripts/cropGenomeByMotifs.py 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30_prefixed.fasta ./0000_DATA/$RegionBoundaries 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30_fixed.fasta
fi
# map output reads to the final consensus sequence with smalt
smalt index -k 15 -s 3 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30_fixed 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30_fixed.fasta
smalt map -x -y 0.7 -j 0 -i 2000 -o 0100_shiver/"${Prefix}"_remap_fixed.sam 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30_fixed ./0020_procreads/$PreprocessedForward ./0020_procreads/$PreprocessedReverse
# sort reads with samtools and index with bamtools
samtools sort 0100_shiver/"${Prefix}"_remap_fixed.sam -o 0100_shiver/"${Prefix}"_remap_fixed.bam
bamtools index -in 0100_shiver/"${Prefix}"_remap_fixed.bam
# call variants with lofreq (not part of shiver), compress with bgzip, index with tabix
docker run --rm -it -v $(pwd)/:/data kboltonlab/lofreq /bin/bash -c "cd ../data; lofreq call --call-indels -f ./0100_shiver/'${Prefix}'_remap_consensus_MinCov_15_30_fixed.fasta -o ./0100_shiver/"${GoldRef%.*}"_shiver.vcf 0100_shiver/'${Prefix}'_remap_fixed.bam"
bgzip -c ./0100_shiver/"${GoldRef%.*}"_shiver.vcf > ./0100_shiver/"${GoldRef%.*}"_shiver.vcf.gz
tabix -fp vcf ./0100_shiver/"${GoldRef%.*}"_shiver.vcf.gz

## Smaltalign
# copy input files into directory
cp ./0000_DATA/$AssemblyRef ./0200_smaltalign # reference sequence
cp ./0020_procreads/"${RawForwardMod%_*.*.*}"_R1_R2.fastq.gz ./0200_smaltalign # interleaved input reads
# enter directory
cd 0200_smaltalign
# source conda.sh, activate conda environment, run Smaltalign, exit conda environment
source /home/hazaihiv/Software/miniconda3/etc/profile.d/conda.sh
conda activate SmaltAlign
/usr/bin/time -q -f "User_time: %U\nSystem_time: %S\nElapsed_time: %E\nMRSS: %M\nFile_outputs: %O\nCPU_percentage: %P\n" -o smaltalign_tm.txt \
    /home/hazaihiv/Software/SmaltAlign/smaltalign_indel.sh -r $AssemblyRef -n 200000
conda deactivate
# create dummy fasta if pipeline fails
if [ ! -f "${RawForwardMod%_*.*.*}"_R1_R2_4_cons.fasta ]
then
    echo "Smaltalign genome doesn't exist. Prepare dummy genome..."
    echo ">smaltalign_error" > "${RawForwardMod%_*.*.*}"_R1_R2_4_cons_fixed.fasta
    printf 'N%.0s' {1..8500} >> "${RawForwardMod%_*.*.*}"_R1_R2_4_cons_fixed.fasta
    echo >> "${RawForwardMod%_*.*.*}"_R1_R2_4_cons_fixed.fasta
else
    echo "Smaltalign genome exists. Trimming..."
    cat "${RawForwardMod%_*.*.*}"_R1_R2_4_cons.fasta | sed 's/:/_/g' > "${RawForwardMod%_*.*.*}"_R1_R2_4_cons_prefixed.fasta
    python3 ~/Analyses/myscripts/cropGenomeByMotifs.py "${RawForwardMod%_*.*.*}"_R1_R2_4_cons_prefixed.fasta ../0000_DATA/$RegionBoundaries "${RawForwardMod%_*.*.*}"_R1_R2_4_cons_fixed.fasta
fi
# remove input and junk files, exit directory
rm $AssemblyRef "${RawForwardMod%_*.*.*}"_R1_R2.fastq.gz *.bam* *.depth #*.vcf
cd ..
# Mapping and variant calling
# map output reads to the final consensus sequence with smalt
smalt index -k 15 -s 3 ./0200_smaltalign/"${RawForwardMod%_*.*.*}"_R1_R2_4_cons_fixed ./0200_smaltalign/"${RawForwardMod%_*.*.*}"_R1_R2_4_cons_fixed.fasta
smalt map -x -y 0.7 -j 0 -i 2000 -o ./0200_smaltalign/"${RawForwardMod%_*.*.*}"_final_mapped.sam ./0200_smaltalign/"${RawForwardMod%_*.*.*}"_R1_R2_4_cons_fixed ./0020_procreads/$PreprocessedForward ./0020_procreads/$PreprocessedReverse
# sort reads with samtools and index with bamtools
samtools sort ./0200_smaltalign/"${RawForwardMod%_*.*.*}"_final_mapped.sam -o ./0200_smaltalign/"${RawForwardMod%_*.*.*}"_final_mapped.bam
bamtools index -in ./0200_smaltalign/"${RawForwardMod%_*.*.*}"_final_mapped.bam
# call variants using lofreq from the final mapping, compress with bgzip, index with tabix
docker run --rm -it -v $(pwd)/:/data kboltonlab/lofreq /bin/bash -c "cd ../data; lofreq call --call-indels -f ./0200_smaltalign/'${RawForwardMod%_*.*.*}'_R1_R2_4_cons_fixed.fasta -o ./0200_smaltalign/'${GoldRef%.*}'_smaltalign.vcf ./0200_smaltalign/'${RawForwardMod%_*.*.*}'_final_mapped.bam"
bgzip -c ./0200_smaltalign/"${GoldRef%.*}"_smaltalign.vcf > ./0200_smaltalign/"${GoldRef%.*}"_smaltalign.vcf.gz
tabix -fp vcf ./0200_smaltalign/"${GoldRef%.*}"_smaltalign.vcf.gz

## viral-ngs
# copy input files into directory
cp ./0000_DATA/$Adapters ./0300_viralngs; cp ./0000_DATA/$AssemblyRef ./0300_viralngs # Adapters and reference sequence
cp ./0000_DATA/$RawForwardMod ./0300_viralngs; cp ./0000_DATA/$RawReverseMod ./0300_viralngs # preprocessed input reads
# enter directory
chmod 777 ./* # temporary
cd 0300_viralngs
# decompress input reads, fix naming to be compatible with later steps
gunzip -d $RawForwardMod; gunzip -d $RawReverseMod
if [[ $Mode == "sim_simref" || $Mode == "sim_SGSref" ]]
then
    cat "${RawForwardMod%.*}" | sed 's/\/1//' > "${RawForwardMod%.*.*}"_r.fastq # this might not be necessary
    cat "${RawReverseMod%.*}" | sed 's/\/2//' > "${RawReverseMod%.*.*}"_r.fastq
else
    cat "${RawForwardMod%.*}" | sed 's/\.1 //' > "${RawForwardMod%.*.*}"_r.fastq # this might need modification
    cat "${RawReverseMod%.*}" | sed 's/\.2 //' > "${RawReverseMod%.*.*}"_r.fastq
fi
# define bam file name
BamReads="${PreprocessedForward%_*_*.*.*}"_R1_R2.bam
# run every step of the viral-ngs pipeline one after another
docker run --rm -it -v $(pwd)/:/data -v /home/hazaihiv/Software/gatk3/:/gatk -v /home/hazaihiv/Software/novoalign/bin/:/novoalign -v /home/hazaihiv/Analyses/myscripts/:/commands quay.io/broadinstitute/viral-ngs /bin/bash -c "apt install time; /usr/bin/time -q -f 'User_time: %U\nSystem_time: %S\nElapsed_time: %E\nMRSS: %M\nFile_outputs: %O\nCPU_percentage: %P\n' -o /data/viralngs_tm.txt bash /commands/vngs.sh $Mode $vngs_prefix $RawForwardMod $RawReverseMod $BamReads $Adapters $AssemblyRef"
# create dummy fasta if pipeline fails
if [ ! -f "${vngs_prefix}"_ordered_imputed_refined.fasta ]
then
    echo "Viralngs genome doesn't exist. Prepare dummy genome..."
    echo ">viralngs_error" > "${vngs_prefix}"_ordered_imputed_refined_fixed.fasta
    printf 'N%.0s' {1..8500} >> "${vngs_prefix}"_ordered_imputed_refined_fixed.fasta
    echo >> "${vngs_prefix}"_ordered_imputed_refined_fixed.fasta
else
    echo "Viralngs genome exists. Trimming..."
    # cut out desired region
    python3 ~/Analyses/myscripts/cropGenomeByMotifs.py "${vngs_prefix}"_ordered_imputed_refined.fasta ../0000_DATA/$RegionBoundaries "${vngs_prefix}"_ordered_imputed_refined_prefixed.fasta
    # change rare ambiguous base characters to N 
    python3 ~/Analyses/myscripts/changeToNbases.py "${vngs_prefix}"_ordered_imputed_refined_prefixed.fasta
    mv genome_mod.fasta "${vngs_prefix}"_ordered_imputed_refined_fixed.fasta
fi
# remove junk and input files, exit directory 
rm $Adapters $AssemblyRef "${RawForwardMod%.*}" "${RawReverseMod%.*}" "${RawForwardMod%.*.*}"_r.fastq "${RawReverseMod%.*.*}"_r.fastq 
cd ..
# map reads to the final genome assembly with smalt, sort with samtools, index with bamtools
smalt index -k 15 -s 3 ./0300_viralngs/"${vngs_prefix}"_ordered_imputed_refined_fixed ./0300_viralngs/"${vngs_prefix}"_ordered_imputed_refined_fixed.fasta
smalt map -x -y 0.7 -j 0 -i 2000 -o ./0300_viralngs/"${PreprocessedForward%_*_*.*.*}"_consmapped.sam ./0300_viralngs/"${vngs_prefix}"_ordered_imputed_refined_fixed ./0020_procreads/$PreprocessedForward ./0020_procreads/$PreprocessedReverse 
samtools sort ./0300_viralngs/"${PreprocessedForward%_*_*.*.*}"_consmapped.sam -o ./0300_viralngs/"${PreprocessedForward%_*_*.*.*}"_consmapped.bam
bamtools index -in ./0300_viralngs/"${PreprocessedForward%_*_*.*.*}"_consmapped.bam
# call variants using lofreq, compress with bgzip, index with tabix
docker run --rm -it -v $(pwd)/:/data kboltonlab/lofreq /bin/bash -c "cd ../data; lofreq call --call-indels -f ./0300_viralngs/'${vngs_prefix}'_ordered_imputed_refined_fixed.fasta -o ./0300_viralngs/'${GoldRef%.*}'_viralngs.vcf ./0300_viralngs/'${PreprocessedForward%_*_*.*.*}'_consmapped.bam"
bgzip -c ./0300_viralngs/"${GoldRef%.*}"_viralngs.vcf > ./0300_viralngs/"${GoldRef%.*}"_viralngs.vcf.gz
tabix -fp vcf ./0300_viralngs/"${GoldRef%.*}"_viralngs.vcf.gz

## V-Pipe
mkdir -p 0400_vpipe/config
# copy data
cp ./0000_DATA/$AssemblyRef ./0400_vpipe/config; cp ./0000_DATA/$Primers ./0400_vpipe/config # reference sequence and primers
cp ./0000_DATA/$VPipeAmpliconProtocol ./0400_vpipe/config # amplicon protocol specification (not working currently)
cp ./0000_DATA/$RawForwardMod ./0400_vpipe; cp ./0000_DATA/$RawReverseMod ./0400_vpipe # original input files (preprocessing done by V-Pipe)
# enter dir
cd 0400_vpipe
# populate 'samples' directory with input reads
mkdir -p ./samples/currsample/currdate/raw_data
mv ./$RawForwardMod ./samples/currsample/currdate/raw_data/reads_R1.fastq.gz
mv ./$RawReverseMod ./samples/currsample/currdate/raw_data/reads_R2.fastq.gz
# initialize project directory
/home/hazaihiv/Software/V-Pipe/V-pipe/init_project.sh 
# override config file with user specification
cp ../0000_DATA/$VPipeConfig .
# # make bed file from primers and reference using my Python script
cd config; python3 ~/Analyses/myscripts/makeMeBed.py $AssemblyRef $Primers; cd ..
# rename reference file to be compatible with the config.yaml
mv ./config/$AssemblyRef ./config/assemblyref.fasta
# create samples.tsv
echo -e "currsample\tcurrdate" > samples.tsv
# calculate then add read length (and amplicon protocol) to samples.tsv
zcat ./samples/currsample/currdate/raw_data/reads_R1.fastq.gz > ./samples/currsample/currdate/raw_data/reads_R1.fastq
Length=$(cat ./samples/currsample/currdate/raw_data/reads_R1.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq | tail -n1)
rm ./samples/currsample/currdate/raw_data/reads_R1.fastq
cat ./samples.tsv | sed "s/$/\t${Length}/" > samples2.tsv
rm samples.tsv; mv samples2.tsv samples.tsv
# fix naming of fastq files
if [[ $Mode != "sim_simref" || $Mode != "sim_SGSref" ]]
then
    zcat ./samples/currsample/currdate/raw_data/reads_R1.fastq.gz | sed 's/\./\//2; s/\./\-/1'  > ./samples/currsample/currdate/raw_data/reads_mod_R1.fastq
    zcat ./samples/currsample/currdate/raw_data/reads_R2.fastq.gz | sed 's/\./\//2; s/\./\-/1'  > ./samples/currsample/currdate/raw_data/reads_mod_R2.fastq
fi
rm ./samples/currsample/currdate/raw_data/reads_R1.fastq.gz ./samples/currsample/currdate/raw_data/reads_R2.fastq.gz
gzip ./samples/currsample/currdate/raw_data/reads_mod_R1.fastq; gzip ./samples/currsample/currdate/raw_data/reads_mod_R2.fastq
# activate conda environment and dry run
source /home/hazaihiv/Software/miniconda3/etc/profile.d/conda.sh
conda activate VpipeMiniconda
./vpipe --dryrun --cores 2
# run V-Pipe
/usr/bin/time -q -f "User_time: %U\nSystem_time: %S\nElapsed_time: %E\nMRSS: %M\nFile_outputs: %O\nCPU_percentage: %P\n" -o vpipe_tm.txt \
    ./vpipe --cores 4
conda deactivate
# save consensus from hmm file (sadly ngshmmalign postprocessing has an issue)
cd results; python3 /home/hazaihiv/Analyses/myscripts/hmmToFasta.py ./currsample/currdate/alignments/currsample-currdate.hmm; cd .. # generates vpipe_consensus_prefixed.fasta
# create dummy fasta if pipeline fails
if [ ! -f results/vpipe_consensus_prefixed.fasta ]
then
    echo "Vpipe genome doesn't exist. Prepare dummy genome..."
    echo ">vpipe_error" > ./results/vpipe_consensus.fasta
    printf 'N%.0s' {1..8500} >> ./results/vpipe_consensus.fasta
    echo >> ./results/vpipe_consensus.fasta
else
    echo "Vpipe genome exists. Trimming..."
    # cut out desired region
    python3 ~/Analyses/myscripts/cropGenomeByMotifs.py ./results/vpipe_consensus_prefixed.fasta ../0000_DATA/$RegionBoundaries ./results/vpipe_consensus.fasta
fi
# remove junk and input files, exit directory 
rm -rf config; cd ..
# map reads to the final genome assembly with smalt, sort with samtools, index with bamtools
smalt index -k 15 -s 3 ./0400_vpipe/results/vpipe_consensus ./0400_vpipe/results/vpipe_consensus.fasta
smalt map -x -y 0.7 -j 0 -i 2000 -o ./0400_vpipe/results/vpipe_final_mapping.sam ./0400_vpipe/results/vpipe_consensus ./0020_procreads/$PreprocessedForward ./0020_procreads/$PreprocessedReverse 
samtools sort ./0400_vpipe/results/vpipe_final_mapping.sam -o ./0400_vpipe/results/vpipe_final_mapping.bam
bamtools index -in ./0400_vpipe/results/vpipe_final_mapping.bam
# call variants using lofreq, compress with bgzip, index with tabix
docker run --rm -it -v $(pwd)/:/data kboltonlab/lofreq /bin/bash -c "cd ../data; lofreq call --call-indels -f ./0400_vpipe/results/vpipe_consensus.fasta -o ./0400_vpipe/results/vpipe_final_mapping.vcf ./0400_vpipe/results/vpipe_final_mapping.bam"
bgzip -c ./0400_vpipe/results/vpipe_final_mapping.vcf > ./0400_vpipe/results/vpipe_final_mapping.vcf.gz
tabix -fp vcf ./0400_vpipe/results/vpipe_final_mapping.vcf.gz

################################
##        Benchmarking       ###
################################

## Populate benchmarking directory
# make dir
mkdir -p 1000_BMresults
# copy preprocessed data
cp 0010_goldenseqs/"$GoldRef" ./1000_BMresults # golden reference sequence
cp 0020_procreads/$PreprocessedForward ./1000_BMresults; cp 0020_procreads/$PreprocessedReverse ./1000_BMresults # preprocessed reads
cp 0030_goldenmap/"${GoldRef%.*}"_goldenmap_sorted_dedupreads.bam ./1000_BMresults # golden mapping
cp 0040_goldenvar/"${GoldRef%.*}".vcf.gz ./1000_BMresults # golden variant file
# copy genome assembly results (mapped reads and consensus genome assembly) and GNU time results
cp  0100_shiver/"${Prefix}"_remap_fixed.bam ./1000_BMresults
cp 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30_fixed.fasta ./1000_BMresults
cp 0100_shiver/"${GoldRef%.*}"_shiver.vcf.gz ./1000_BMresults
cp 0200_smaltalign/"${RawForwardMod%_*.*.*}"_final_mapped.bam ./1000_BMresults
cp 0200_smaltalign/"${GoldRef%.*}"_smaltalign.vcf.gz ./1000_BMresults
cp 0200_smaltalign/"${RawForwardMod%_*.*.*}"_R1_R2_4_cons_fixed.fasta ./1000_BMresults
cp 0300_viralngs/"${PreprocessedForward%_*_*.*.*}"_consmapped.bam ./1000_BMresults
cp 0300_viralngs/"${vngs_prefix}"_ordered_imputed_refined_fixed.fasta ./1000_BMresults
cp 0300_viralngs/"${GoldRef%.*}"_viralngs.vcf.gz ./1000_BMresults
cp 0400_vpipe/results/vpipe_final_mapping.bam ./1000_BMresults
cp 0400_vpipe/results/vpipe_consensus.fasta ./1000_BMresults
cp 0400_vpipe/results/vpipe_final_mapping.vcf.gz ./1000_BMresults
cp 0100_shiver/shiver_tm.txt ./1000_BMresults
cp 0200_smaltalign/smaltalign_tm.txt ./1000_BMresults
cp 0300_viralngs/viralngs_tm.txt ./1000_BMresults
cp 0400_vpipe/vpipe_tm.txt ./1000_BMresults
# Rename files
cd 1000_BMresults
mv  "${Prefix}"_remap_fixed.bam shiver_reads.bam
mv "${Prefix}"_remap_consensus_MinCov_15_30_fixed.fasta shiver_genome.fasta
mv "${GoldRef%.*}"_shiver.vcf.gz shiver_variants.vcf.gz
mv "${RawForwardMod%_*.*.*}"_final_mapped.bam smaltalign_reads.bam
mv "${RawForwardMod%_*.*.*}"_R1_R2_4_cons_fixed.fasta smaltalign_genome.fasta
mv "${GoldRef%.*}"_smaltalign.vcf.gz smaltalign_variants.vcf.gz
mv "${PreprocessedForward%_*_*.*.*}"_consmapped.bam viralngs_reads.bam
mv "${vngs_prefix}"_ordered_imputed_refined_fixed.fasta viralngs_genome.fasta
mv "${GoldRef%.*}"_viralngs.vcf.gz viralngs_variants.vcf.gz
mv vpipe_final_mapping.bam vpipe_reads.bam
mv vpipe_consensus.fasta vpipe_genome.fasta
mv vpipe_final_mapping.vcf.gz vpipe_variants.vcf.gz

## Fix files
# fix naming of bam files
samtools view -h shiver_reads.bam | sed 's+\.1\t+\/1\t+g' | sed 's+\.2\t+\/2\t+g' | samtools view -bS > shiver_reads_m.bam
samtools view -h smaltalign_reads.bam | sed 's+\.1\t+\/1\t+g' | sed 's+\.2\t+\/2\t+g' | samtools view -bS > smaltalign_reads_m.bam
samtools view -h viralngs_reads.bam | sed 's+\.1\t+\/1\t+g' | sed 's+\.2\t+\/2\t+g' | samtools view -bS > viralngs_reads_m.bam
samtools view -h vpipe_reads.bam | sed 's+\.1\t+\/1\t+g' | sed 's+\.2\t+\/2\t+g' | samtools view -bS > vpipe_reads_m.bam
samtools view -h "${GoldRef%.*}"_goldenmap_sorted_dedupreads.bam | sed 's+\.1\t+\/1\t+g' | sed 's+\.2\t+\/2\t+g' | samtools view -bS > "${GoldRef%.*}"_goldenmap_sorted_dedupreads_m.bam
# decompress read files, fix naming of fastq files
gzip -d $PreprocessedForward; gzip -d $PreprocessedReverse
cat "${PreprocessedForward%.*.*}".fastq | sed 's+\.1+\/1+g' | sed 's+\.2+\/2+g' > "${PreprocessedForward%.*.*}"_m.fastq
cat "${PreprocessedReverse%.*.*}".fastq | sed 's+\.1+\/1+g' | sed 's+\.2+\/2+g' > "${PreprocessedReverse%.*.*}"_m.fastq

## Variant analysis
cat $GoldRef <(echo) shiver_genome.fasta smaltalign_genome.fasta viralngs_genome.fasta vpipe_genome.fasta  > "assemblies.fasta"
mafft --localpair --maxiterate 1000 assemblies.fasta > assemblies_aligned.fasta
gzip -d "${GoldRef%.*}".vcf.gz; gzip -d shiver_variants.vcf.gz; gzip -d smaltalign_variants.vcf.gz; gzip -d viralngs_variants.vcf.gz; gzip -d vpipe_variants.vcf.gz
python /home/hazaihiv/Analyses/myscripts/variantAnalysis.py assemblies_aligned.fasta  "${GoldRef%.*}".vcf shiver_variants.vcf smaltalign_variants.vcf viralngs_variants.vcf vpipe_variants.vcf

## Mapping analysis
# get name of reads
echo "Start mapping analysis:"
samtools view  "${GoldRef%.*}"_goldenmap_sorted_dedupreads_m.bam | cut -f1  | sort | uniq > names.txt
# get file with the sequence coordinates
echo "Here comes that I am a fool, but I know that pretty well..."
for i in "${GoldRef%.*}"_goldenmap_sorted_dedupreads_m.bam shiver_reads_m.bam smaltalign_reads_m.bam viralngs_reads_m.bam vpipe_reads_m.bam; do java -jar ~/Software/jvarkit/jvarkit/dist/jvarkit.jar samgrep -f names.txt $i > "${i%.*}".sam; done
java -jar ~/Software/jvarkit/jvarkit/dist/cmpbams.jar -F "${GoldRef%.*}"_goldenmap_sorted_dedupreads_m.sam shiver_reads_m.sam smaltalign_reads_m.sam viralngs_reads_m.sam vpipe_reads_m.sam > compbams.txt
# Identify mapping errors:
python3 ~/Analyses/myscripts/compBamLiftover.py assemblies_aligned.fasta compbams.txt 
# tidy up
#rm *.sam names.txt compbams.txt

## QUAST report
# Run QUAST
docker run --rm -it -v $(pwd)/:/data staphb/quast /bin/bash -c "quast.py shiver_genome.fasta smaltalign_genome.fasta viralngs_genome.fasta vpipe_genome.fasta --reference '$GoldRef' --pe1 '${PreprocessedForward%.*.*}'_m.fastq --pe2 '${PreprocessedReverse%.*.*}'_m.fastq --bam shiver_reads_m.bam,smaltalign_reads_m.bam,viralngs_reads_m.bam,vpipe_reads_m.bam --ref-bam '${GoldRef%.*}'_goldenmap_sorted_dedupreads_m.bam"

## get some basic results from the simulations
if [[ $Mode == "sim_simref" ]]
then
    # generations
    echo -n "Generations: " | tee bm_stats.txt
    echo "$Generation" | tee -a bm_stats.txt
    # length of benchmarkign reference
    echo -n "Length of benchmarking reference: " | tee -a bm_stats.txt
    awk '/^>/ {seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' $GoldRef | tee -a bm_stats.txt
    # distance from HXB2
    echo -n "Distance from HXB2: " | tee -a bm_stats.txt
    cat ../0000_DATA/NCBI_HIV1_HXB2.fasta $GoldRef > refs.fasta
    python3 ~/Analyses/myscripts/alignMAFFT.py refs.fasta
    python3 ~/Analyses/myscripts/multifastaDiversity.py multifasta_aligned.fasta | tee -a bm_stats.txt
    rm diversity_matrix.csv fasta_upper.fasta multifasta_aligned.fasta refs.fasta
    # QS diversity
    echo -n "Average pairwise diversity of quasispecies: " | tee -a bm_stats.txt
    python3 ~/Analyses/myscripts/cropGenomeByMotifs.py ../0000_DATA/$SeqSet ../0000_DATA/HXB2_genome_boundaries.fasta QS_fixed.fasta
    python3 ~/Analyses/myscripts/alignMAFFT.py QS_fixed.fasta
    python3 ~/Analyses/myscripts/multifastaDiversity.py multifasta_aligned.fasta | tee -a bm_stats.txt
    rm diversity_matrix.csv fasta_upper.fasta multifasta_aligned.fasta QS_fixed.fasta
fi

## clean up temporary files
cd ..
rm -rf 0010_goldenseqs 0020_procreads 0030_goldenmap 0040_goldenvar 0100_shiver 0200_smaltalign 0300_viralngs 0400_vpipe