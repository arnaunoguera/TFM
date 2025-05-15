#!/bin/bash

## WARNING
# Requires conda env seqtk

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -S <step> -P <passwd> ]"
    echo
    echo "Options:"
    echo "  -i  PLASMID_DIR   Input dir (with redundant plasmids fasta file)"
    echo "              Default: BASE_DIR/PLASMIDS"
    echo "  -o  MOBMESS_DIR  Output directory to store results. Default: BASE_DIR/MOBMESS"
    echo "  -h          Display this help message"
}

#Pots haver de fer abans: export TMPDIR=BASE_DIR/MOBMESS/tmp/

# Get options to run different parts of the strainphlan pipeline
while getopts i:o:t:P:s:h flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) results_dir=${OPTARG};;
        t) threads=${OPTARG};;
        P) passwd=${OPTARG};;
        s) step=${OPTARG};;
        h) show_help
           exit 0;;
        *) show_help
           exit 1;;
    esac
done


###### DEFINE VARIABLES TO RUN the program ######


if [ -z "$threads" ]; then
  threads=1
fi

if [ -z "$input_dir" ]; then
  input_dir="BASE_DIR/PLASMIDS"
fi

if [ -z "$results_dir" ]; then
  results_dir="BASE_DIR/MOBMESS"
fi

if [ -z "$step" ]; then
  step=1
fi


# Main script


echo "Filtering clean plasmids from the redundant plasmids file"
# Get the fasta file with the clean plasmids
seqtk subseq "$input_dir/predicted_plasmids_redundant.fa" "$input_dir/clean_plasmid_list.txt" > "$input_dir/clean_plasmids_redundant.fa"

echo "Filtering clean plasmids from the redundant plasmids info file"
# Filter the plasmids info file
# First save the header
head -n 1 $input_dir/redundant_plasmid_final_info.txt > $input_dir/clean_redundant_plasmids_info.txt
grep -f $input_dir/clean_plasmid_list.txt $input_dir/redundant_plasmid_final_info.txt >> $input_dir/clean_redundant_plasmids_info.txt

echo "Filtering clean plasmids from the redundant plasmids alignment file"
# Filter the plasmid alignment file from mmseq2
grep -v -F -f $input_dir/non_plasmid_list.txt $results_dir/predicted_plasmids_pairwise_alignments.txt > $results_dir/clean_redundant_plasmids_pairwise_alignments.txt