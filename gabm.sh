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
    echo "m     Mode ('emp_ref' - empirical data with known sequences,'emp_noref' - empirical data with unknown sequences,'sim' - simulated data)"
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
        m) # Enter a name
            Mode=$OPTARG
            Accepted=("emp_ref" "emp_noref" "sim")
            if [[ " ${Accepted[*]} " =~ " ${Mode} " ]]
            then
                echo "Mode: $Mode"
            else
                echo "Error: Invalid mode"
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
if [[ $Mode == "emp_ref" || $Mode == "emp_noref" ]]
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
            ForwardInput=$ForwardInput\.gz # rename bash variable
            ReverseInput=$ReverseInput\.gz
        else
            echo "There were no .fastq files in the directory!"
            echo "Put input files into the ./0000_DATA directory according to the manual. Exiting..."
            exit 0
        fi
    fi
fi


############################
## Benchmarking baseline  ##
############################

## QS sequence generation with santa-sim
if [[ $Mode == "sim" ]]
then
    cd 0000_DATA; java -jar /home/hazaihiv/Software/santa-sim/dist/santa.jar $SantaConf; mv alignment_1.fasta "${Master%.*}_QSsim.fasta"; cd ..
    SeqSet="${Master%.*}_QSsim.fasta" # define sequence set for the 'sim' mode
fi

## Processing QS sequences
mkdir -p 0010_goldenseqs # make dir
if [[ $Mode != "emp_noref" ]]
then
    # multiple alignment
    mafft --localpair --maxiterate 1000 ./0000_DATA/"${SeqSet}" > ./0010_goldenseqs/"${SeqSet%.*}"_aligned.fasta
    # mafft ./0000_DATA/"${SeqSet}" > ./0010_goldenseqs/"${SeqSet%.*}"_aligned.fasta # just for development purposes
    # consensus
    sudo chmod +777 -R ../* # solve permission issue with output files
    # call consensus from the QS sequence set
    docker run --rm -it -v $(pwd)/:/data biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1 /bin/bash -c "em_cons -sequence ./0010_goldenseqs/'${SeqSet%.*}'_aligned.fasta -outseq ./0010_goldenseqs/'${SeqSet%.*}'_aligned_consensus.fasta"
    samtools faidx ./0010_goldenseqs/"${SeqSet%.*}"_aligned_consensus.fasta # index consensus
fi

## Define golden reference sequence
if [[ $Mode != "emp_noref" ]]
then
    GoldRef="${SeqSet%.*}"_aligned_consensus.fasta # consensus of known sequences
else
    sudo chmod +777 -R ../*
    cp ./0000_DATA/$GoldRef ./0010_goldenseqs # arbitrary sequence (results might indicate outlier assemblies between genome assembly software)
fi

## Simulation of PE NGS reads
if [[ $Mode == "sim" ]]
then
    # generate sequences with ART
	docker run --rm -it -v $(pwd)/:/data 19da3c67bd0edccc378d135e79b7ce8a376d3aabccfb087c97cb340920ff1eed /bin/bash -c "cd 0000_DATA; art_illumina -sam -i '$SeqSet' -l 250 -s 1 -f '$SimCov' -ss MS -o '${SeqSet%.*}' -m 700 -nf 0 -p -qL 20 -qU 40 -rs 1680244882"
    # remove junk files
    sudo rm ./0000_DATA/"${SeqSet%_*.*}"_QSsim1.aln ./0000_DATA/"${SeqSet%_*.*}"_QSsim2.aln ./0000_DATA/"${SeqSet%.*}".sam
    # move results to the data directory
    mv ./0000_DATA/"${SeqSet%_*.*}"_QSsim1.fq ./0000_DATA/"${SeqSet%_*.*}"_QSsim_R1.fastq; mv ./0000_DATA/"${SeqSet%_*.*}"_QSsim2.fq ./0000_DATA/"${SeqSet%_*.*}"_QSsim_R2.fastq
    # compress the simulated reads
    gzip ./0000_DATA/"${SeqSet%_*.*}"_QSsim_R1.fastq; gzip ./0000_DATA/"${SeqSet%_*.*}"_QSsim_R2.fastq
    # name read input variables
    ForwardInput="${SeqSet%_*.*}"_QSsim_R1.fastq.gz
    ReverseInput="${SeqSet%_*.*}"_QSsim_R2.fastq.gz
fi

## Preprocessing PE NGS reads 
# make dir
mkdir -p 0020_procreads
# QC
fastqc -o 0020_procreads ./0000_DATA/$ForwardInput
fastqc -o 0020_procreads ./0000_DATA/$ReverseInput
# Read trimming and filtering with trimmomatic
TrimmomaticPE ./0000_DATA/$ForwardInput ./0000_DATA/$ReverseInput ./0020_procreads/"${ForwardInput%.*.*}"_paired.fastq.gz ./0020_procreads/"${ForwardInput%.*.*}"_unpaired.fastq.gz ./0020_procreads/"${ReverseInput%.*.*}"_paired.fastq.gz ./0020_procreads/"${ReverseInput%.*.*}"_unpaired.fastq.gz ILLUMINACLIP:./0000_DATA/$Adapters:2:10:7:1:true MINLEN:50 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20
gzip ./0020_procreads/"${ForwardInput%.*.*}"_paired.fastq; gzip ./0020_procreads/"${ReverseInput%.*.*}"_paired.fastq
# name preprocessed reads
ForwardReads="${ForwardInput%.*.*}"_paired.fastq.gz
ReverseReads="${ReverseInput%.*.*}"_paired.fastq.gz
# interleave forward and reverse reads (for smaltalign):
reformat.sh in1=./0000_DATA/$ForwardInput in2=./0000_DATA/$ReverseInput out=./0020_procreads/"${ForwardInput%_*.*.*}"_R1_R2.fastq.gz

## Golden mapping
# make dir
mkdir -p 0030_goldenmap
########################################################################################
# # ngshmmalign (writing of processed reads is not working) 
# gzip -d ./0020_procreads/${ForwardReads%.*}; gzip -d ./0020_procreads/${ReverseReads%.*}
# docker run --rm -it -v $(pwd)/:/data quay.io/biocontainers/ngshmmalign:0.1.1--ha04c180_4 /bin/bash -c "ngshmmalign -h;cd data;ngshmmalign -v --input ./0020_procreads/${ForwardReads%.*} ./0020_procreads/${ReverseReads%.*} -r ./0030_goldenmap/"${SeqSet%.*}"_aligned.fasta -o ./0030_goldenmap/"${SeqSet%.*}"_aligned.sam"
# gzip ./0020_procreads/$ForwardReads; gzip ./0020_procreads/$ForwardReads
##########################################################################################
# index consensus and map empirical reads with smalt to the consensus
smalt index -k 15 -s 3 ./0010_goldenseqs/"${GoldRef%.*}" ./0010_goldenseqs/$GoldRef
smalt map -x -y 0.7 -j 0 -i 2000 -o ./0030_goldenmap/"${GoldRef%.*}"_goldenmap.sam ./0010_goldenseqs/"${GoldRef%.*}" ./0020_procreads/$ForwardReads ./0020_procreads/$ReverseReads 
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



#################################
###      Genome assembly      ###
#################################

# make directories
mkdir -p 0100_shiver
mkdir -p 0200_smaltalign
mkdir -p 0300_viralngs
mkdir -p 0400_vpipe/config

## Shiver
# copy input files to shiver directory
cp ./0000_DATA/$RefAlignment ./0100_shiver # alignment of reference sequences
cp ./0000_DATA/$Adapters ./0100_shiver; cp ./0000_DATA/$Primers ./0100_shiver # adapters and primers
cp ./0000_DATA/$ShiverConfig ./0100_shiver; cp ./0000_DATA/pipeline.conf ./0100_shiver # shiver config files
cp ./0000_DATA/$ForwardInput ./0100_shiver; cp ./0000_DATA/$ReverseInput ./0100_shiver # processed input reads (trimmomatic inside the shiver pipeline, but same specification)
# go inside directory
cd 0100_shiver
# run shiver pipeline
/usr/bin/time -v -o shiver_tm.txt \
    docker run --rm -it -v $(pwd)/:/data fazekasda/shiver:22.04 full
# remove input files and exit directory
rm $RefAlignment $Adapters $ShiverConfig pipeline.conf $Primers $ForwardInput $ReverseInput
cd ..
# shiver consensus file contains 2 sequences, remove one without reference imputation
cat 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30.fasta | awk '/>'"${Prefix}"'_remap_consensus/ {getline; while(!/>/) {getline}} 1' > 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30_fixed.fasta
# call variants with lofreq (not part of shiver), compress with bgzip, index with tabix
docker run --rm -it -v $(pwd)/:/data kboltonlab/lofreq /bin/bash -c "cd ../data; lofreq call --call-indels -f ./0100_shiver/'${Prefix}'_remap_consensus_MinCov_15_30_fixed.fasta -o ./0100_shiver/'${GoldRef%.*}'_shiver.vcf ./0100_shiver/'${Prefix}'_remap.bam"
bgzip -c ./0100_shiver/"${GoldRef%.*}"_shiver.vcf > ./0100_shiver/"${GoldRef%.*}"_shiver.vcf.gz
tabix -fp vcf ./0100_shiver/"${GoldRef%.*}"_shiver.vcf.gz

## Smaltalign
# copy input files into directory
cp ./0000_DATA/$AssemblyRef ./0200_smaltalign # reference sequence
cp ./0020_procreads/"${ForwardInput%_*.*.*}"_R1_R2.fastq.gz ./0200_smaltalign # interleaved input reads
# enter directory
cd 0200_smaltalign
# source conda.sh, activate conda environment, run Smaltalign, exit conda environment
source /home/hazaihiv/Software/miniconda3/etc/profile.d/conda.sh
conda activate SmaltAlign
/usr/bin/time -v -o smaltalign_tm.txt \
    /home/hazaihiv/Software/SmaltAlign/smaltalign.sh -r $AssemblyRef
conda deactivate
# remove input and junk files, exit directory
rm $AssemblyRef "${ForwardInput%_*.*.*}"_R1_R2.fastq.gz *.bam* *.depth #*.vcf
cd ..
# Mapping and variant calling
if [[ $Mode != "emp_noref" ]]
then
    # map output reads to the final consensus sequence with smalt
    smalt index -k 15 -s 3 ./0200_smaltalign/"${ForwardInput%_*.*.*}"_R1_R2_4_cons ./0200_smaltalign/"${ForwardInput%_*.*.*}"_R1_R2_4_cons.fasta
    smalt map -x -y 0.7 -j 0 -i 2000 -o ./0200_smaltalign/"${ForwardInput%_*.*.*}"_final_mapped.sam ./0200_smaltalign/"${ForwardInput%_*.*.*}"_R1_R2_4_cons ./0020_procreads/$ForwardReads ./0020_procreads/$ReverseReads
    # sort reads with samtools and index with bamtools
    samtools sort ./0200_smaltalign/"${ForwardInput%_*.*.*}"_final_mapped.sam -o ./0200_smaltalign/"${ForwardInput%_*.*.*}"_final_mapped.bam
    bamtools index -in ./0200_smaltalign/"${ForwardInput%_*.*.*}"_final_mapped.bam
    # call variants using lofreq from the final mapping, compress with bgzip, index with tabix
    docker run --rm -it -v $(pwd)/:/data kboltonlab/lofreq /bin/bash -c "cd ../data; lofreq call --call-indels -f ./0200_smaltalign/'${ForwardInput%_*.*.*}'_R1_R2_4_cons.fasta -o ./0200_smaltalign/'${GoldRef%.*}'_smaltalign.vcf ./0200_smaltalign/'${ForwardInput%_*.*.*}'_final_mapped.bam"
    bgzip -c ./0200_smaltalign/"${GoldRef%.*}"_smaltalign.vcf > ./0200_smaltalign/"${GoldRef%.*}"_smaltalign.vcf.gz
    tabix -fp vcf ./0200_smaltalign/"${GoldRef%.*}"_smaltalign.vcf.gz
else
    # map output reads to the final consensus sequence with smalt
    smalt index -k 15 -s 3 ./0200_smaltalign/"${ForwardInput%_*_*.*.*}"_4_cons ./0200_smaltalign/"${ForwardInput%_*_*.*.*}"_4_cons.fasta
    smalt map -x -y 0.7 -j 0 -i 2000 -o ./0200_smaltalign/"${ForwardInput%_*.*.*}"_final_mapped.sam ./0200_smaltalign/"${ForwardInput%_*_*.*.*}"_4_cons ./0020_procreads/$ForwardReads ./0020_procreads/$ReverseReads
    # sort reads with samtools and index with bamtools
    samtools sort ./0200_smaltalign/"${ForwardInput%_*.*.*}"_final_mapped.sam -o ./0200_smaltalign/"${ForwardInput%_*.*.*}"_final_mapped.bam
    bamtools index -in ./0200_smaltalign/"${ForwardInput%_*.*.*}"_final_mapped.bam
    # call variants using lofreq from the final mapping, compress with bgzip, index with tabix
    docker run --rm -it -v $(pwd)/:/data kboltonlab/lofreq /bin/bash -c "cd ../data; lofreq call --call-indels -f ./0200_smaltalign/'${ForwardInput%_*_*.*.*}'_4_cons.fasta -o ./0200_smaltalign/'${GoldRef%.*}'_smaltalign.vcf ./0200_smaltalign/'${ForwardInput%_*.*.*}'_final_mapped.bam"
    bgzip -c ./0200_smaltalign/"${GoldRef%.*}"_smaltalign.vcf > ./0200_smaltalign/"${GoldRef%.*}"_smaltalign.vcf.gz
    tabix -fp vcf ./0200_smaltalign/"${GoldRef%.*}"_smaltalign.vcf.gz
fi

## viral-ngs
# copy input files into directory
cp ./0000_DATA/$Adapters ./0300_viralngs; cp ./0000_DATA/$AssemblyRef ./0300_viralngs # Adapters and reference sequence
cp ./0000_DATA/$ForwardInput ./0300_viralngs; cp ./0000_DATA/$ReverseInput ./0300_viralngs # preprocessed input reads
# enter directory
chmod 777 ./* # temporary
cd 0300_viralngs
# decompress input reads, fix naming to be compatible with later steps
gunzip -d $ForwardInput; gunzip -d $ReverseInput
if [[ $Mode == "sim" ]]
then
    cat "${ForwardInput%.*}" | sed 's/\/1//' > "${ForwardInput%.*.*}"_r.fastq # this might not be necessary
    cat "${ReverseInput%.*}" | sed 's/\/2//' > "${ReverseInput%.*.*}"_r.fastq
else
    cat "${ForwardInput%.*}" | sed 's/\.1 //' > "${ForwardInput%.*.*}"_r.fastq # this might need modification
    cat "${ReverseInput%.*}" | sed 's/\.2 //' > "${ReverseInput%.*.*}"_r.fastq
fi
# define bam file name
BamReads="${ForwardReads%_*_*.*.*}"_R1_R2.bam
# run every step of the viral-ngs pipeline one after another
if [[ $Mode != "sim" ]]
then
    /usr/bin/time -v -o viralngs_tm.txt \
        docker run --rm -it -v $(pwd)/:/data -v /home/hazaihiv/Software/gatk3/:/gatk -v /home/hazaihiv/Software/novoalign/bin/:/novoalign quay.io/broadinstitute/viral-ngs /bin/bash -c "cd /data; /opt/viral-ngs/source/read_utils.py fastq_to_bam --sampleName $vngs_prefix '${ForwardInput%.*.*}'_r.fastq '${ReverseInput%.*.*}'_r.fastq init.bam; /opt/viral-ngs/source/read_utils.py sort_bam init.bam $BamReads coordinate; rm init.bam; /opt/viral-ngs/source/assembly.py assemble_spades $BamReads $Adapters:2:10:7:1:true '${vngs_prefix}'.fasta; /opt/viral-ngs/source/assembly.py order_and_orient '${vngs_prefix}'.fasta $AssemblyRef '${vngs_prefix}'_ordered.fasta; /opt/viral-ngs/source/assembly.py impute_from_reference '${vngs_prefix}'_ordered.fasta $AssemblyRef '${vngs_prefix}'_ordered_imputed.fasta --index; /opt/viral-ngs/source/assembly.py refine_assembly '${vngs_prefix}'_ordered_imputed.fasta $BamReads '${vngs_prefix}'_ordered_imputed_refined.fasta"
else
    /usr/bin/time -v -o viralngs_tm.txt \
        docker run --rm -it -v $(pwd)/:/data -v /home/hazaihiv/Software/gatk3/:/gatk -v /home/hazaihiv/Software/novoalign/bin/:/novoalign quay.io/broadinstitute/viral-ngs /bin/bash -c "cd /data; /opt/viral-ngs/source/read_utils.py fastq_to_bam --sampleName $vngs_prefix '${ForwardInput%.*.*}'_r.fastq '${ReverseInput%.*.*}'_r.fastq init.bam; /opt/viral-ngs/source/read_utils.py sort_bam init.bam $BamReads coordinate; rm init.bam; /opt/viral-ngs/source/assembly.py assemble_trinity $BamReads $Adapters:2:10:7:1:true '${vngs_prefix}'.fasta; /opt/viral-ngs/source/assembly.py order_and_orient '${vngs_prefix}'.fasta $AssemblyRef '${vngs_prefix}'_ordered.fasta; /opt/viral-ngs/source/assembly.py impute_from_reference '${vngs_prefix}'_ordered.fasta $AssemblyRef '${vngs_prefix}'_ordered_imputed.fasta --index; /opt/viral-ngs/source/assembly.py refine_assembly '${vngs_prefix}'_ordered_imputed.fasta $BamReads '${vngs_prefix}'_ordered_imputed_refined.fasta"
fi
remove junk and input files, exit directory 
rm $Adapters $AssemblyRef "${ForwardInput%.*}" "${ReverseInput%.*}" "${ForwardInput%.*.*}"_r.fastq "${ReverseInput%.*.*}"_r.fastq 
cd ..
# map reads to the final genome assembly with smalt, sort with samtools, index with bamtools
smalt index -k 15 -s 3 ./0300_viralngs/"${vngs_prefix}"_ordered_imputed_refined ./0300_viralngs/"${vngs_prefix}"_ordered_imputed_refined.fasta
smalt map -x -y 0.7 -j 0 -i 2000 -o ./0300_viralngs/"${ForwardReads%_*_*.*.*}"_consmapped.sam ./0300_viralngs/"${vngs_prefix}"_ordered_imputed_refined ./0020_procreads/$ForwardReads ./0020_procreads/$ReverseReads 
samtools sort ./0300_viralngs/"${ForwardReads%_*_*.*.*}"_consmapped.sam -o ./0300_viralngs/"${ForwardReads%_*_*.*.*}"_consmapped.bam
bamtools index -in ./0300_viralngs/"${ForwardReads%_*_*.*.*}"_consmapped.bam
# call variants using lofreq, compress with bgzip, index with tabix
docker run --rm -it -v $(pwd)/:/data kboltonlab/lofreq /bin/bash -c "cd ../data; lofreq call --call-indels -f ./0300_viralngs/'${vngs_prefix}'_ordered_imputed_refined.fasta -o ./0300_viralngs/'${GoldRef%.*}'_viralngs.vcf ./0300_viralngs/'${ForwardReads%_*_*.*.*}'_consmapped.bam"
bgzip -c ./0300_viralngs/"${GoldRef%.*}"_viralngs.vcf > ./0300_viralngs/"${GoldRef%.*}"_viralngs.vcf.gz
tabix -fp vcf ./0300_viralngs/"${GoldRef%.*}"_viralngs.vcf.gz

## V-Pipe
mkdir -p 0400_vpipe/config
# copy data
cp ./0000_DATA/$AssemblyRef ./0400_vpipe/config; cp ./0000_DATA/$Primers ./0400_vpipe/config # reference sequence and primers
cp ./0000_DATA/$VPipeAmpliconProtocol ./0400_vpipe/config # amplicon protocol specification (not working currently)
cp ./0000_DATA/$ForwardInput ./0400_vpipe; cp ./0000_DATA/$ReverseInput ./0400_vpipe # original input files (preprocessing done by V-Pipe)
# enter dir
cd 0400_vpipe
# populate 'samples' directory with input reads
mkdir -p ./samples/currsample/currdate/raw_data
mv ./$ForwardInput ./samples/currsample/currdate/raw_data/reads_R1.fastq.gz
mv ./$ReverseInput ./samples/currsample/currdate/raw_data/reads_R2.fastq.gz
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
if [[ $Mode != "sim" ]]
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
/usr/bin/time -v -o vpipe_tm.txt \
    ./vpipe --cores 4
conda deactivate
save consensus from hmm file (sadly ngshmmalign postprocessing has an issue)
cd results; python3 /home/hazaihiv/Analyses/myscripts/hmmToFasta.py ./currsample/currdate/alignments/currsample-currdate.hmm; cd .. # generates vpipe_consensus.fasta
# remove junk and input files, exit directory 
rm -rf config; cd ..
# map reads to the final genome assembly with smalt, sort with samtools, index with bamtools
smalt index -k 15 -s 3 ./0400_vpipe/results/vpipe_consensus ./0400_vpipe/results/vpipe_consensus.fasta
smalt map -x -y 0.7 -j 0 -i 2000 -o ./0400_vpipe/results/vpipe_final_mapping.sam ./0400_vpipe/results/vpipe_consensus ./0020_procreads/$ForwardReads ./0020_procreads/$ReverseReads 
samtools sort ./0400_vpipe/results/vpipe_final_mapping.sam -o ./0400_vpipe/results/vpipe_final_mapping.bam
bamtools index -in ./0400_vpipe/results/vpipe_final_mapping.bam
# call variants using lofreq, compress with bgzip, index with tabix
docker run --rm -it -v $(pwd)/:/data kboltonlab/lofreq /bin/bash -c "cd ../data; lofreq call --call-indels -f ./0400_vpipe/results/vpipe_consensus.fasta -o ./0400_vpipe/results/vpipe_final_mapping.vcf ./0400_vpipe/results/vpipe_final_mapping.bam"
bgzip -c ./0400_vpipe/results/vpipe_final_mapping.vcf > ./0400_vpipe/results/vpipe_final_mapping.vcf.gz
tabix -fp vcf ./0400_vpipe/results/vpipe_final_mapping.vcf.gz

################################
##        Benchmarking       ###
################################

# GoldRef="NCBI_HIV1_HXB2_QSsim_aligned_consensus.fasta"
# ForwardInput="NCBI_HIV1_HXB2_QSsim_R1.fastq.gz"
# ReverseInput="NCBI_HIV1_HXB2_QSsim_R2.fastq.gz"
# ForwardReads="NCBI_HIV1_HXB2_QSsim_R1_paired.fastq.gz"
# ReverseReads="NCBI_HIV1_HXB2_QSsim_R2_paired.fastq.gz"

# # only for testing purposes
# ForwardInput="5VM_R1.fastq.gz"
# ReverseInput="5VM_R2.fastq.gz"
# ForwardReads="${ForwardInput%.*.*}"_paired.fastq.gz
# ReverseReads="${ReverseInput%.*.*}"_paired.fastq.gz
# GoldRef="${SeqSet%.*}"_aligned_consensus.fasta # consensus of known sequences

## Populate benchmarking directory
# make dir
mkdir -p 1000_BMresults
# copy preprocessed data
cp 0010_goldenseqs/"$GoldRef" ./1000_BMresults # golden reference sequence
cp 0020_procreads/$ForwardReads ./1000_BMresults; cp 0020_procreads/$ReverseReads ./1000_BMresults # preprocessed reads
cp 0030_goldenmap/"${GoldRef%.*}"_goldenmap_sorted_dedupreads.bam ./1000_BMresults # golden mapping
cp 0040_goldenvar/"${GoldRef%.*}".vcf.gz ./1000_BMresults # golden variant file
# copy genome assembly results (mapped reads and consensus genome assembly)
cp 0100_shiver/"${Prefix}"_remap.bam ./1000_BMresults
cp 0100_shiver/"${Prefix}"_remap_consensus_MinCov_15_30_fixed.fasta ./1000_BMresults
cp 0100_shiver/"${GoldRef%.*}"_shiver.vcf.gz ./1000_BMresults
cp 0200_smaltalign/"${ForwardInput%_*.*.*}"_final_mapped.bam ./1000_BMresults
cp 0200_smaltalign/"${GoldRef%.*}"_smaltalign.vcf.gz ./1000_BMresults
if [[ $Mode != "emp_noref" ]]
then
    cp 0200_smaltalign/"${ForwardInput%_*.*.*}"_R1_R2_4_cons.fasta ./1000_BMresults
else
    cp 0200_smaltalign/"${ForwardInput%_*_*.*.*}"_4_cons.fasta ./1000_BMresults
fi
cp 0300_viralngs/"${ForwardReads%_*_*.*.*}"_consmapped.bam ./1000_BMresults
cp 0300_viralngs/"${vngs_prefix}"_ordered_imputed_refined.fasta ./1000_BMresults
cp 0300_viralngs/"${GoldRef%.*}"_viralngs.vcf.gz ./1000_BMresults
cp 0400_vpipe/results/vpipe_final_mapping.bam ./1000_BMresults
cp 0400_vpipe/results/vpipe_consensus.fasta ./1000_BMresults
cp 0400_vpipe/results/vpipe_final_mapping.vcf.gz ./1000_BMresults
# Rename files
cd 1000_BMresults
mv "${Prefix}"_remap.bam shiver_reads.bam
mv "${Prefix}"_remap_consensus_MinCov_15_30_fixed.fasta shiver_genome.fasta
mv "${GoldRef%.*}"_shiver.vcf.gz shiver_variants.vcf.gz
mv "${ForwardInput%_*.*.*}"_final_mapped.bam smaltalign_reads.bam
if [[ $Mode != "emp_noref" ]]
then
    mv "${ForwardInput%_*.*.*}"_R1_R2_4_cons.fasta smaltalign_genome.fasta
else
    mv "${ForwardInput%_*_*.*.*}"_4_cons.fasta smaltalign_genome.fasta
fi
mv "${GoldRef%.*}"_smaltalign.vcf.gz smaltalign_variants.vcf.gz
mv "${ForwardReads%_*_*.*.*}"_consmapped.bam viralngs_reads.bam
mv "${vngs_prefix}"_ordered_imputed_refined.fasta viralngs_genome.fasta
mv "${GoldRef%.*}"_viralngs.vcf.gz viralngs_variants.vcf.gz
mv vpipe_final_mapping.bam vpipe_reads.bam
mv vpipe_consensus.fasta vpipe_genome.fasta
mv vpipe_final_mapping.vcf.gz vpipe_variants.vcf.gz

## Variant analysis
cat $GoldRef <(echo) shiver_genome.fasta smaltalign_genome.fasta viralngs_genome.fasta vpipe_genome.fasta > "assemblies.fasta"
mafft --localpair --maxiterate 1000 assemblies.fasta > assemblies_aligned.fasta
gzip -d "${GoldRef%.*}".vcf.gz; gzip -d shiver_variants.vcf.gz; gzip -d smaltalign_variants.vcf.gz; gzip -d viralngs_variants.vcf.gz; gzip -d vpipe_variants.vcf.gz
python /home/hazaihiv/Analyses/myscripts/variantAnalysis.py assemblies_aligned.fasta  "${GoldRef%.*}".vcf shiver_variants.vcf smaltalign_variants.vcf viralngs_variants.vcf vpipe_variants.vcf

## QUAST report
# fix naming of bam files
samtools view -h shiver_reads.bam | sed 's+\.1\t+\/1\t+g' | sed 's+\.2\t+\/2\t+g' | samtools view -bS > shiver_reads_m.bam
samtools view -h smaltalign_reads.bam | sed 's+\.1\t+\/1\t+g' | sed 's+\.2\t+\/2\t+g' | samtools view -bS > smaltalign_reads_m.bam
samtools view -h viralngs_reads.bam | sed 's+\.1\t+\/1\t+g' | sed 's+\.2\t+\/2\t+g' | samtools view -bS > viralngs_reads_m.bam
samtools view -h vpipe_reads.bam | sed 's+\.1\t+\/1\t+g' | sed 's+\.2\t+\/2\t+g' | samtools view -bS > vpipe_reads_m.bam
samtools view -h "${GoldRef%.*}"_goldenmap_sorted_dedupreads.bam | sed 's+\.1\t+\/1\t+g' | sed 's+\.2\t+\/2\t+g' | samtools view -bS > "${GoldRef%.*}"_goldenmap_sorted_dedupreads_m.bam
# decompress read files, fix naming of fastq files
gzip -d $ForwardReads; gzip -d $ReverseReads
cat "${ForwardReads%.*.*}".fastq | sed 's+\.1+\/1+g' | sed 's+\.2+\/2+g' > "${ForwardReads%.*.*}"_m.fastq
cat "${ReverseReads%.*.*}".fastq | sed 's+\.1+\/1+g' | sed 's+\.2+\/2+g' > "${ReverseReads%.*.*}"_m.fastq
# Run QUAST
docker run --rm -it -v $(pwd)/:/data staphb/quast /bin/bash -c "quast.py shiver_genome.fasta smaltalign_genome.fasta viralngs_genome.fasta vpipe_genome.fasta --reference '$GoldRef' --pe1 '${ForwardReads%.*.*}'_m.fastq --pe2 '${ReverseReads%.*.*}'_m.fastq --bam shiver_reads_m.bam --bam smaltalign_reads_m.bam --bam viralngs_reads_m.bam --bam vpipe_reads_m.bam --ref-bam '${GoldRef%.*}'_goldenmap_sorted_dedupreads.bam"
# clean up, exit directory
rm  "$GoldRef" "${ForwardReads%.*.*}".fastq "${ReverseReads%.*.*}".fastq "${ForwardReads%.*.*}"_m.fastq "${ReverseReads%.*.*}"_m.fastq "${GoldRef%.*}"_goldenmap_sorted_dedupreads.bam "${GoldRef%.*}"_goldenmap_sorted_dedupreads_m.bam "${GoldRef%.*}".vcf
cd ..
