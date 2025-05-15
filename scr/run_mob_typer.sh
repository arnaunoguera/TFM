#!/bin/bash


# Script to: run mob-suite to predict the mobility and incompatibility groups of the plasmids
#######WARNING######
# This sctipt neds the conda environment mob_suite activated

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -S <step> -P <passwd> ]"
    echo
    echo "Options:"
    echo "  -p  PLASMIDS_DOIR   Input dir (with redundant plasmids fasta file)"
    echo "                      Default: BASE_DIR/PLASMIDS"
    echo "  -o  MOB-SUITE DIR   Directory for results of mob-suite. Default: BASE_DIR/MOB_SUITE"
    echo "  -s  STEP            Step to run. Default: 1. Options:"
    echo "                      - 1: Split the plasmid fasta file into n parts to run one after the other (not to overcharge memory) (needs conda env seqkit)"  
    echo "                      - 2: Run mob-suite to predict the mobility and incompatibility groups of the plasmids for each part"
    echo "                      - 3: Merge all the information into a single file"
    echo "  -t  THREADS         Number of threads to use (default: 1)"
    echo "  -n  N_PARTITIONS   Number of partitions to split the input file (default: 30)"
    echo "  -h          Display this help message"
}

# Get options to run different parts of the strainphlan pipeline
while getopts p:P:v:b:m:s:t:n:h flag
do
    case "${flag}" in
        p) plasmid_dir=${OPTARG};;
        o) MOB_dir=${OPTARG};;
        s) step=${OPTARG};;
        t) threads=${OPTARG};;
        n) n_partitions=${OPTARG};;
        h) show_help
           exit 0;;
        *) show_help
           exit 1;;
    esac
done


###### DEFINE VARIABLES TO RUN PLASX ######


if [ -z "$plasmid_dir" ]; then
  plasmid_dir="BASE_DIR/PLASMIDS"
fi

if [ -z "$MOB_dir" ]; then
  MOB_dir="BASE_DIR/MOB_SUITE"
fi

if [ -z "$step" ]; then
  step="1"
fi

if [ -z "$threads" ]; then
  threads=1
fi

if [ -z "$n_partitions" ]; then
  n_partitions=30
fi

### MAIN SCRIPT

if [ $step == "1" ]; then
    echo "Splitting the plasmid fasta file into $n_partitions parts"
    seqkit split2 -p $n_partitions -O $plasmid_dir/parallel_runs/ -j $threads $plasmid_dir/predicted_plasmids_redundant.fa
fi

if [ $step == "2" ]; then
    #run mob-typer for each split file
    for i in $(seq -f "%03g" 1 $n_partitions); do
        echo "Running mob-typer for partition $i"
        # split files are in parallel_runs/predicted_plasmids_redundant.part_019.fa
        mob_typer -d /mnt/synology/DATABASES/MOB_SUITE/ -o $MOB_dir/mobtyper_results_$i.txt -i $plasmid_dir/parallel_runs/predicted_plasmids_redundant.part_$i.fa -n $threads --multi
    done
fi

if [ $step == "3" ]; then
    echo "Merging all the mob-typer results into a single file"
    # print the header of the first file
    head -n 1 $MOB_dir/mobtyper_results_001.txt > $MOB_dir/mobtyper_results.txt
    # concat the rest of the files without the header
    tail -q -n +2 $MOB_dir/mobtyper_results_*.txt >> $MOB_dir/mobtyper_results.txt

    # remove temp files 
    #rm $MOB_dir/mobtyper_results_*.txt
fi