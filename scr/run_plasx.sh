#!/bin/bash

#######WARNING######
# This sctipt neds the conda environment plasx activated
# The COG and Pfam annotations are done in these specific version because they're the ones required for PlasX

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -S <step> -P <passwd> ]"
    echo
    echo "Options:"
    echo "  -i  INPUT   Input directory (directory with the anvi'o results, with one folder per sample inside)"
    echo "              Default: BASE_DIR/ANVIO/POP"
    echo "  -o  OUTPUT  Output directory to store results. Default: BASE_DIR/PLASX/POP"
    echo "  -t  THREADS Number of threads (default: 1)"
    echo "  -s  SAMPLE LIST "
    echo "              File with a list of the samples to run this script for. "
    echo "              Default: BASE_DIR/SAMPLE_LISTS/POP.txt"
    echo "  -P  PSSWD   Password to run sudo commands"
    echo "  -h          Display this help message"
}

# Get options to run different parts of the strainphlan pipeline
while getopts i:o:t:s:P:S:h flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) results_dir=${OPTARG};;
        t) threads=${OPTARG};;
        s) sample_list=${OPTARG};;
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
  input_dir="BASE_DIR/ANVIO/POP"
fi

if [ -z "$results_dir" ]; then
  results_dir="BASE_DIR/PLASX/POP"
fi

if [ -z "$sample_list" ]; then
  sample_list="BASE_DIR/SAMPLE_LISTS/POP.txt"
fi


#####MAIN SCRIPT######

if [ -d $input_dir ]; then
  echo "Input directory is $input_dir"
else
  echo "ERROR: Couldn't find input directory. Exitting..."
  exit 1
fi

time for sample in $(cat $sample_list);
do
  

  echo 
  echo "Starting process for sample $sample ..."
  if ! [ -d "$results_dir/${sample}/" ]; then
    mkdir -pv "$results_dir/${sample}"
  else
   echo "$results_dir/${sample}/ already exists, skipping to the next sample..."
   continue
  fi
  

  echo
  # Annotate genes using the database of de novo gene families
  echo "Annotating genes using the database of de novo gene families..."
  # - This command requires a high amount of RAM (at least ~60Gb). If your machine has low RAM, then you need to set the `--splits` flag to a high number.
  #   This will split the de novo families into chunks, reducing the max RAM usage. E.g. if you have only ~8Gb, we recommend setting `--splits` to 32 or higher.
  plasx search_de_novo_families \
      -g "$input_dir/$sample/$sample-gene-calls.txt" \
      -o "$results_dir/$sample/$sample-de-novo-families.txt" \
      --threads $threads \
      --overwrite

  # Use PlasX to classify contigs as plasmid or non-plasmid sequences, based on the annotations from previous steps
  echo
  echo "Predicting plasmids with PlasX..."
  plasx predict \
    -a "$input_dir/$sample/$sample-cogs-and-pfams.txt" "$results_dir/$sample/$sample-de-novo-families.txt" \
    -g "$input_dir/$sample/$sample-gene-calls.txt" \
    -o "$results_dir/$sample/$sample-scores.txt" \
    --overwrite

  # Clear swap
  echo 'metalab00' | sudo -S swapoff -a; sudo -S swapon -a

done
