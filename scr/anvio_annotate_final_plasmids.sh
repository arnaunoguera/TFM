#!/bin/bash

#######WARNING######
# This sctipt neds the conda environment anvio-8 activated
# The COG and Pfam annotations are done in these specific version because they're the ones required for PlasX

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -S <step> -P <passwd> ]"
    echo
    echo "Options:"
    echo "  -i  INPUT   Input directory with the final plasmids fasta file"
    echo "              Default: BASE_DIR/PLASMIDS"
    echo "  -o  OUTPUT  Output directory to store results. Default: BASE_DIR/ANVIO/FINAL_PLASMIDS"
    echo "  -t  THREADS Number of threads (default: 1)"
    echo "  -S  STEP    Step to run. Options: [1, 2, 3, 4, 5, 6, 7, all]"
    echo "              1: Create an anvio contigs database from the fasta file"
    echo "              2: Export gene calls (including amino acid sequences) to text file"
    echo "              3: Annotate COGs"
    echo "              4: Annotate Pfams"
    echo "              5: Export functions to text file"
    echo "              6: Run HMMs"
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
pfam_dir="/mnt/synology/DATABASES/Pfam_v37.2"
COG_dir="/mnt/synology/DATABASES/COG_2020"


if [ -z "$threads" ]; then
  threads=1
fi

if [ -z "$input_dir" ]; then
  input_dir="BASE_DIR/PLASMIDS"
fi

if [ -z "$results_dir" ]; then
  results_dir="BASE_DIR/ANVIO/FINAL_PLASMIDS"
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


log_file="$results_dir/anvio.log"

  
if [ "$step" = "1" ] || [ "$step" = "all" ]; then
echo
# Create an anvio contigs database from the fasta file
echo "Creating an anvio contigs database from contigs.fasta..."
# - The `-L 0` parameter ensures that contigs remain intact and aren't split
anvi-gen-contigs-database -L 0 -T $threads --project-name nonredundant_plasmids -f "$input_dir/nonredundant_plasmids.fa" \
    -o "$results_dir/nonredundant_plasmids.db" 2>&1 | tee -a "$log_file"
fi

if [ "$step" = "2" ] || [ "$step" = "all" ]; then
echo
# Export gene calls (including amino acid sequences) to text file
echo "Exporting gene calls to text nonredundant_plasmids-gene-calls.txt..."
anvi-export-gene-calls --gene-caller prodigal -c "$results_dir/nonredundant_plasmids.db" \
    -o "$results_dir/nonredundant_plasmids-gene-calls.txt" 2>&1 | tee -a "$log_file"
fi

if [ "$step" = "3" ] || [ "$step" = "all" ]; then
echo
# Annotate COGs
echo "Annotationg COGs..."
anvi-run-ncbi-cogs -T $threads --cog-version COG20 --cog-data-dir $COG_dir -c "$results_dir/nonredundant_plasmids.db" 2>&1 | tee -a "$log_file"
fi

if [ "$step" = "4" ] || [ "$step" = "all" ]; then
echo
# Annotate Pfams
echo "Annotationg Pfam..."
anvi-run-pfams -T $threads --pfam-data-dir $pfam_dir -c "$results_dir/nonredundant_plasmids.db" 2>&1 | tee -a "$log_file"
fi

if [ "$step" = "5" ] || [ "$step" = "all" ]; then
echo
# Export functions to text file
echo "Exporting functions to text file..."
anvi-export-functions --annotation-sources "COG20_FUNCTION,Pfam" -c "$results_dir/nonredundant_plasmids.db" \
    -o "$results_dir/nonredundant_plasmids-cogs-and-pfams.txt" 2>&1 | tee -a "$log_file"
fi

if [ "$step" = "6" ] || [ "$step" = "all" ]; then
echo 
# Running HMMs
echo "Running HMMs..."
anvi-run-hmms -c "$results_dir/nonredundant_plasmids.db" --num-threads "$threads" 2>&1 | tee -a "$log_file"
fi
