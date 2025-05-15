#!/bin/bash

#######WARNING######
# This sctipt neds the conda environment rgi activated

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -d <database> [ -S <step> -P <passwd> ]"
    echo
    echo "Options:"
    echo "  -i  INPUT   Input directory (directory with the plasmid database fasta file)"
    echo "              Default: BASE_DIR/PLASMIDS"
    echo "  -o  OUTPUT  Output directory to store results. Default: BASE_DIR/CARD"
    echo "  -t  THREADS Number of threads (default: 1)"
    echo "  -h          Display this help message"
}

# Get options to run different parts of the strainphlan pipeline
while getopts i:o:t:d:S:h flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) results_dir=${OPTARG};;
        t) threads=${OPTARG};;
        P) passwd=${OPTARG};;
        h) show_help
           exit 0;;
        *) show_help
           exit 1;;
    esac
done


###### DEFINE VARIABLES TO RUN PLASX ######


if [ -z "$threads" ]; then
  threads=1
fi

if [ -z "$input_dir" ]; then
  input_dir="BASE_DIR/PLASMIDS"
fi

if [ -z "$results_dir" ]; then
  results_dir="BASE_DIR/CARD"
fi

#####MAIN SCRIPT######
 
rgi main -i $input_dir/nonredundant_plasmids.fa -o $results_dir/CARD_hits -t contig -n $threads -d plasmid --split_prodigal_jobs --clean