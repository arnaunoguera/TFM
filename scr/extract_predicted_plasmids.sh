#!/bin/bash


# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -a <anvio_dir> -s <sample_list> [ -p <preset> -c ]"
    echo
    echo "Options:"
    echo "  -i  INPUT   Input directory (directory with the results of PLASX, with one folder per sample inside)"
    echo "              Default: BASE_DIR/PLASX/POP"
    echo "  -o  OUTPUT  Output directory to store assemblies. Default: BASE_DIR/PLASMIDS/POP"
    echo "  -a  Directory with the anvio name-corrected fasta files (default: BASE_DIR/ANVIO/POP)"
    echo "  -s  SAMPLE LIST "
    echo "              File with a list of the samples to run this script for. "
    echo "              Default: BASE_DIR/samplelist.txt"
    echo "  -S  STEP    Step to run. Default: 1. Options:"
    echo "                  - 1: Extract predicted plasmids"
    echo "                  - 2: For each sample, extract the anvio annotations for the predicted plasmids"
    echo "  -h          Display this help message"    
}

# Get options to run different parts of the strainphlan pipeline
while getopts i:o:a:s:S:h flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) results_dir=${OPTARG};;
        s) sample_list=${OPTARG};;
        a) anvio_dir=${OPTARG};;
        S) step=${OPTARG};;
        h) show_help
           exit 0;;
        *) show_help
           exit 1;;
    esac
done


if [ -z "$input_dir" ]; then
  input_dir="BASE_DIR/PLASX/POP"
fi

if [ -z "$results_dir" ]; then
  results_dir="BASE_DIR/PLASMIDS/POP"
fi

if [ -z "$sample_list" ]; then
  sample_list="BASE_DIR/SAMPLE_LISTS/POP.txt"
fi

if [ -z "$anvio_dir" ]; then
  anvio_dir="BASE_DIR/ANVIO/POP.txt"
fi

if [ -z "$step" ]; then
  step="1"
fi

#####MAIN SCRIPT######
if [ -d $input_dir ]; then
  echo "Input directory is $input_dir"
else
  echo "ERROR: Couldn't find input directory. Exitting..."
  exit 1
fi

for sample in $(cat $sample_list);
do

  if [ "$step" == "1" ]; then
    echo "Extracting predicted plasmids for sample $sample"
    if ! [ -d "$results_dir/${sample}/" ]; then
      mkdir -pv "$results_dir/${sample}"
    else
    echo "$results_dir/${sample}/ already exists, skipping to the next sample..."
    continue
    fi

    # Extract predicted plasmids (where the value > 0.5)
    # settinc LC numeric to C so awk interprets them as numbers properly
    LC_NUMERIC=C awk 'NR > 1 {if ($2 > 0.5) print $0}' "$input_dir/$sample/$sample-scores.txt" > "$results_dir/$sample/$sample-predicted-plasmids.txt"

    # Get a list with just the names to extract them 
    awk '{print $1}' "$results_dir/$sample/$sample-predicted-plasmids.txt" > "$input_dir/$sample/$sample-predicted-plasmids.txt"
    
    #Extract the fasta sequences from the list of predicted plasmids
    seqtk subseq "$anvio_dir/$sample/$sample.fixed_contigs.fa" "$input_dir/$sample/$sample-predicted-plasmids.txt" > "$results_dir/$sample/$sample-predicted-plasmids.fa"
  fi 

  if [ "$step" == "2" ]; then
    echo "Extracting annotations for sample $sample"

    # Extract the gene calls for the predicted plasmids and get the gene called id list  
    grep -f "$input_dir/$sample/$sample-predicted-plasmids.txt" "$anvio_dir/$sample/$sample-gene-calls.txt" | cut -f 1 > "$results_dir/$sample/$sample-gene-callers_id.txt"
    # Then save the selected gene calls but adding the sample anme to the gene call id (to be able to merge them later)
    head -n 1 "$anvio_dir/$sample/$sample-gene-calls.txt" > "$results_dir/$sample/$sample-gene-calls.txt"
    grep -f "$input_dir/$sample/$sample-predicted-plasmids.txt" "$anvio_dir/$sample/$sample-gene-calls.txt" | awk -v sample="$sample" '{print sample"_"$0}' >> "$results_dir/$sample/$sample-gene-calls.txt"

    # Extract the cogs and pfams for the predicted plasmids, using the gene_callers_id file
    head -n 1 "$anvio_dir/$sample/$sample-cogs-and-pfams.txt" > "$results_dir/$sample/$sample-cogs-and-pfams.txt"
    awk -v sample="$sample" 'FNR==NR { ids[$1]; next } $1 in ids {print sample"_"$0}' "$results_dir/$sample/$sample-gene-callers_id.txt" "$anvio_dir/$sample/$sample-cogs-and-pfams.txt" >> "$results_dir/$sample/$sample-cogs-and-pfams.txt"

  fi 
done
