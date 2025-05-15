#!/bin/bash

#######WARNING######
# This sctipt neds the conda environment anvio8 activated
# The COG and Pfam annotations are done in these specific version because they're the ones required for PlasX

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -S <step> -P <passwd> ]"
    echo
    echo "Options:"
    echo "  -i  INPUT   Input directory (directory with the assembled metagenomes, with one folder per sample inside)"
    echo "              Default: BASE_DIR/MEGAHIT/POP"
    echo "  -o  OUTPUT  Output directory to store results. Default: BASE_DIR/ANVIO/POP"
    echo "  -t  THREADS Number of threads (default: 1)"
    echo "  -s  SAMPLE LIST "
    echo "              File with a list of the samples to run this script for. "
    echo "              Default: BASE_DIR/SAMPLE_LISTS/POP.txt"
    echo "  -P  PSSWD   Password to run sudo commands"
    echo "  -S  STEP    Step to run. Options: [1, 2, 3, 4, 5, 6, 7, all]"
    echo "              1: Reformat FASTA file to simplify names"
    echo "              2: Create an anvio contigs database from the fasta file"
    echo "              3: Export gene calls (including amino acid sequences) to text file"
    echo "              4: Annotate COGs"
    echo "              5: Annotate Pfams"
    echo "              6: Export functions to text file"
    echo "              7: Run HMMs"
    echo "              all: Run all steps (default)"
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
        S) step=${OPTARG};;
        h) show_help
           exit 0;;
        *) show_help
           exit 1;;
    esac
done

###### DEFINE VARIABLES TO RUN ANVI'O ######
pfam_dir="/mnt/synology/DATABASES/Pfam_v32"
COG_dir="/mnt/synology/DATABASES/COG_2014"


if [ -z "$threads" ]; then
  threads=1
fi

if [ -z "$input_dir" ]; then
  input_dir="BASE_DIR/MEGAHIT/POP"
fi

if [ -z "$results_dir" ]; then
  results_dir="BASE_DIR/ANVIO/POP"
fi

if [ -z "$sample_list" ]; then
  sample_list="BASE_DIR/SAMPLE_LISTS/POP.txt"
fi

if [ -z "$step" ]; then
  step="all"
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

  log_file="$results_dir/$sample/$sample.log"

  if [ "$step" = "1" ] || [ "$step" = "all" ]; then
    echo
    # Reformat FASTA file to simplify names and keep only contigs longer than 1kbp if it wasn't done before
    echo "Reformatting FASTA file to simplify contig names..."
    anvi-script-reformat-fasta "$input_dir/$sample/$sample.contigs.fa" \
      -o "$results_dir/$sample/$sample.fixed_contigs.fa" -l 1000 --simplify-names --prefix $sample \
      --report "$results_dir/$sample/$sample.name_conversions.txt" 2>&1 | tee -a "$log_file"
  fi

  if [ "$step" = "2" ] || [ "$step" = "all" ]; then
    echo
    # Create an anvio contigs database from the fasta file
    echo "Creating an anvio contigs database from contigs.fasta..."
    # - The `-L 0` parameter ensures that contigs remain intact and aren't split
    anvi-gen-contigs-database -L 0 -T $threads --project-name $sample -f "$results_dir/$sample/$sample.fixed_contigs.fa" \
      -o "$results_dir/$sample/$sample.db" 2>&1 | tee -a "$log_file"
  fi

  if [ "$step" = "3" ] || [ "$step" = "all" ]; then
    echo
    # Export gene calls (including amino acid sequences) to text file
    echo "Exporting gene calls to text $sample-gene-calls.txt..."
    anvi-export-gene-calls --gene-caller prodigal -c "$results_dir/$sample/$sample.db" \
      -o "$results_dir/$sample/$sample-gene-calls.txt" 2>&1 | tee -a "$log_file"
  fi

  if [ "$step" = "4" ] || [ "$step" = "all" ]; then
    echo
    # Annotate COGs
    echo "Annotationg COGs..."
    anvi-run-ncbi-cogs -T $threads --cog-version COG14 --cog-data-dir $COG_dir -c "$results_dir/$sample/$sample.db" 2>&1 | tee -a "$log_file"
  fi

  if [ "$step" = "5" ] || [ "$step" = "all" ]; then
    echo
    # Annotate Pfams
    echo "Annotationg Pfam..."
    anvi-run-pfams -T $threads --pfam-data-dir $pfam_dir -c "$results_dir/$sample/$sample.db" 2>&1 | tee -a "$log_file"
  fi

  if [ "$step" = "6" ] || [ "$step" = "all" ]; then
    echo
    # Export functions to text file
    echo "Exporting functions to text file..."
    anvi-export-functions --annotation-sources "COG14_FUNCTION,Pfam" -c "$results_dir/$sample/$sample.db" \
      -o "$results_dir/$sample/$sample-cogs-and-pfams.txt" 2>&1 | tee -a "$log_file"
  fi

  if [ "$step" = "7" ] || [ "$step" = "all" ]; then
    echo 
    # Running HMMs
    echo "Running HMMs..."
    anvi-run-hmms -c "$results_dir/$sample/$sample.db" --num-threads "$threads" 2>&1 | tee -a "$log_file"
  fi

  # Clear swap
  echo $passwd | sudo -S swapoff -a; sudo -S swapon -a

done
