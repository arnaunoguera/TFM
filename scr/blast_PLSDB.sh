#!/bin/bash


# Script to: run blastn alignments of each plasmid against the PLSDB database to check if certain contigs are plasmids for sure
#######WARNING######
# This sctipt neds the conda environment blast activated

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -S <step> -P <passwd> ]"
    echo
    echo "Options:"
    echo "  -i  INPUT   Input dir (with redundant plasmids fasta file)"
    echo "              Default: BASE_DIR/PLASMIDS"
    echo "  -o  OUTPUT  Output directory to store the pairwise alignment results. Default: BASE_DIR/PLSDB_BLAST"
    echo "  -t  THREADS Number of threads (default: 1)"
    # echo "  -s  STEP    Step to run. Default: 1. Options:"
    # echo "                  - 1: Run the pairwise alignments with mmseq2"
    # echo "                  - 2: Split input file into 20 parts to parallelize"
    # echo "                  - 3: Annotate circularit based on ViralVerify and pairwise alignments results"
    # echo "                  - 4: Merge parallelized results"
    echo "  -d DATABASE PLSDB blast database (default: /mnt/synology/DATABASES/PLSDB/PLSDB_2024_05_31_v2.fasta)"
    echo "  -h          Display this help message"
}

# Get options to run different parts of the strainphlan pipeline
while getopts i:o:t:v:n:d:s:m:r:h flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) results_dir=${OPTARG};;
        t) threads=${OPTARG};;
        s) step=${OPTARG};;
        d) PLSDB=${OPTARG};;
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
  results_dir="BASE_DIR/PLSDB_BLAST"
fi

if [ -z "$step" ]; then
  step="1"
fi

if [ -z "$PLSDB" ]; then 
  PLSDB="/mnt/synology/DATABASES/PLSDB/PLSDB_2024_05_31_v2.fasta"
fi


#### MAIN SCRIPT 

# Run blast with all the database. Parameters:
# -evalue 1e-8: only keep hits with very significant e-value
# -outfmt "6 qaccver saccver pident qcovus length qlen qstart qend sstart send evalue bitscore": custom tab-separated format.
# Special interest in qcovus (query coverage per subject): the percentage of the query sequence that matches the subject sequence (adding up all local alignments)
# pident: percentage of identity
# num_alignments 5: number of hits to keep per query sequence
# perc_identity 90: only keep hits with at least 90% identity
# num_threads: number of threads to use

echo "Running blastn against the PLSDB database..."

blastn -query $input_dir/predicted_plasmids_redundant.fa -db $PLSDB -evalue 1e-8 -outfmt "6 qaccver saccver pident qcovus length qlen qstart qend sstart send evalue bitscore" -num_alignments 5 \
    -perc_identity 90 -num_threads $threads > $results_dir/blast_results.tsv

## Note that stitle is the description

echo "Done!"
echo "Extracting hits with at least 90% coverage..."

# Get the list of plasmids and their hits in the database if the coverage of the query (qcovus) is enough (90%)
LC_NUMERIC=C awk '$4 >= 90 {print $0}' $results_dir/blast_results.tsv > $results_dir/blast_results_COV90.tsv
LC_NUMERIC=C awk '$4 >= 90 {print $1 "\t" $2}' $results_dir/blast_results.tsv | sort | uniq > $results_dir/PLSDB_hits.tsv


# Get the list of plasmids that have hits in the database
cut -f 1 $results_dir/PLSDB_hits.tsv | sort | uniq > $results_dir/PLSDB_hits_plasmids.txt

