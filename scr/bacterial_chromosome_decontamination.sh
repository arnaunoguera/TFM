#!/bin/bash


# Script to: detect chromosomal contigs in the plasmid fasta file by aligning against database of reference bacterial genomes (with filtering of plasmid sequences)
#######WARNING######
# This sctipt neds the conda environment ncbi_datasets activated

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -S <step> -P <passwd> ]"
    echo
    echo "Options:"
    echo "  -i  INPUT   Input dir (with redundant plasmids fasta file)"
    echo "              Default: BASE_DIR/PLASMIDS"
    echo "  -o  OUTPUT  Output directory to store results. Default: BASE_DIR/BACTERIAL_CHROMOSOME_FILTERING"
    echo "  -t  THREADS Number of threads (default: 1)"
    echo "  -s  STEP    Step to run. Default: 1. Options:"
    echo "                  - 1: Download reference bacterial genomes from NCBI (requires conda environment ncbi_datasets)"
    echo "                  - 2: Build chromosomal database by removing plasmid sequences (requires conda environment seqtk)"
    echo "                  - 3: Build blast database with the reference chromosomes (requires conda environment blast)"
    echo "                  - 4: Run blast to detect chromosomal contigs in the plasmid fasta file (requires conda environment blast)"
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
  results_dir="BASE_DIR/BACTERIAL_CHROMOSOME_FILTERING"
fi

if [ -z "$step" ]; then
  step="1"
fi

### MAIN SCRIPT ###

if [ $step == "1" ]; then
    echo "Downloading reference bacterial genomes from NCBI"

    # Download the dehydrated reference genomes
    datasets download genome taxon 2,archaea --assembly-level chromosome,complete --assembly-version latest --exclude-atypical --reference --dehydrated --chromosomes 1,chromosome \
        --filename $results_dir/BACTERIAL_REFERENCE_GENOMES/ncbi_genomes_dehydrated.zip 

    # Unzip the reference genomes
    unzip $results_dir/BACTERIAL_REFERENCE_GENOMES/ncbi_genomes_dehydrated.zip -d $results_dir/BACTERIAL_REFERENCE_GENOMES/ncbi_genomes_dehydrated

    # Rehydrate the reference genomes
    datasets rehydrate --directory $results_dir/BACTERIAL_REFERENCE_GENOMES/ncbi_genomes_dehydrated --max-workers 30

    mkdir -p $results_dir/BACTERIAL_REFERENCE_GENOMES/FASTA

    # extract fasta files into a single folder
    mv  $results_dir/BACTERIAL_REFERENCE_GENOMES/ncbi_genomes_dehydrated/ncbi_dataset/data/*/*.fna $results_dir/BACTERIAL_REFERENCE_GENOMES/FASTA/

    # remove ncbi directories
    rm $results_dir/BACTERIAL_REFERENCE_GENOMES/ncbi_genomes_dehydrated.zip 
    rm -r $results_dir/BACTERIAL_REFERENCE_GENOMES/ncbi_genomes_dehydrated

fi 


if [ $step == "2" ]; then

    echo "Building downloaded database of reference bacterial genomes"
    # Put all the files in a single fasta file
    cat $results_dir/BACTERIAL_REFERENCE_GENOMES/FASTA/* > $results_dir/BACTERIAL_REFERENCE_GENOMES/all_genomes.fa

    echo "Extracting chromosome sequences from the reference genomes"
    # Extract the first sequence header from each file, which is the bacterial genome in all cases
    # I checked and there are no first-line headers with the word plasmid in the header
    head -n 1 -q $results_dir/BACTERIAL_REFERENCE_GENOMES/all_genomes.fa > $results_dir/BACTERIAL_REFERENCE_GENOMES/first_headers.txt

    # Also extract those that contain the word chromosome
    grep -ih "chromosome" $results_dir/BACTERIAL_REFERENCE_GENOMES/all_genomes.fa > $results_dir/BACTERIAL_REFERENCE_GENOMES/chromosome_headers.txt

    # Join both files 
    cat $results_dir/BACTERIAL_REFERENCE_GENOMES/first_headers.txt $results_dir/BACTERIAL_REFERENCE_GENOMES/chromosome_headers.txt | grep "[^>]*" -io | sort | uniq > $results_dir/BACTERIAL_REFERENCE_GENOMES/genome_headers.txt

    # remove intermediates
    rm $results_dir/BACTERIAL_REFERENCE_GENOMES/first_headers.txt $results_dir/BACTERIAL_REFERENCE_GENOMES/chromosome_headers.txt

    echo "Filtering plasmid sequences from the reference genomes"
    # Filter the fasta file to keep those that are not plasmids
    seqtk subseq $results_dir/BACTERIAL_REFERENCE_GENOMES/all_genomes.fa $results_dir/BACTERIAL_REFERENCE_GENOMES/genome_headers.txt > $results_dir/BACTERIAL_REFERENCE_GENOMES/reference_chromosomes.fa

    rm $results_dir/BACTERIAL_REFERENCE_GENOMES/all_genomes.fa

fi

if [ "$step" == "3" ]; then

    echo "Building blast database with the reference chromosomes"
    # Build a blast database with the reference chromosomes
    makeblastdb -in $results_dir/BACTERIAL_REFERENCE_GENOMES/reference_chromosomes.fa -dbtype nucl

fi


if [ "$step" == "4" ]; then

    echo "Running blast to detect chromosomal contigs in the plasmid fasta file"
    # Run blast with the plasmid sequences against the reference chromosomes
    blastn -query $input_dir/predicted_plasmids_redundant.fa -db $results_dir/BACTERIAL_REFERENCE_GENOMES/reference_chromosomes.fa -evalue 1e-8 -num_alignments 5 \
        -outfmt "6 qaccver saccver pident qcovus length qlen qstart qend sstart send evalue bitscore" -num_threads $threads -perc_identity 90 > $results_dir/blast_results.tsv

    # Get the list of plasmids and their hits in the database if the coverage of the query (qcovus) is enough (90%)
    LC_NUMERIC=C awk '$4 >= 90 {print $0}' $results_dir/blast_results.tsv > $results_dir/blast_results_COV90.tsv
    LC_NUMERIC=C awk '$4 >= 90 {print $1 "\t" $2}' $results_dir/blast_results.tsv | sort | uniq > $results_dir/bacterial_chromosome_hits.tsv


    # Get the list of plasmids that have hits in the database
    cut -f 1 $results_dir/bacterial_chromosome_hits.tsv | sort | uniq > $results_dir/bacterial_chromosome_contigs.txt


fi