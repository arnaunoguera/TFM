#!/bin/bash

#######WARNING######
# This sctipt neds the conda environment mobmess activated

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -S <step> -P <passwd> ]"
    echo
    echo "Options:"
    echo "  -i  INPUT   Input dir (with redundant plasmids fasta file)"
    echo "              Default: BASE_DIR/PLASMIDS"
    echo "  -o  OUTPUT  Output directory to store results. Default: BASE_DIR/MOBMESS"
    echo "  -t  THREADS Number of threads (default: 1)"
    echo "  -P  PSSWD   Password to run sudo commands"
    echo "  -s  STEP    Step to run. Default: all. Options:"
    echo "                  - 1: Create file with circularity information for all plasmids"
    # echo "                  - 2: Run mmseq2 to get ANI file before running mobmess (requires conda env mmseq2)"
    echo "                  - 2: Run mobmess systems with precoomputed ANI file (requires conda env mobmess)"
    echo "                  - 3: Get deduplicated list of plasmids"
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

if [ -d $input_dir ]; then
  echo "Input directory is $input_dir"
else
  echo "ERROR: Couldn't find input directory. Exitting..."
  exit 1
fi

if [ ! -d $results_dir ]; then
  mkdir -p $results_dir
fi

if [ "$step" == "1" ]; then

    # Create file annotating circularity for all plasmids
    echo "Creating file with circularity information for all plasmids..."

    awk 'NR > 1 {print $1 "\t" ($6 == "circular" ? "1" : "0")}' $input_dir/clean_redundant_plasmids_info.txt > $results_dir/plasmids_circular.txt

fi

# if [ "$step" == "2" ]; then

#     # This commented chunk was to run alignment with nucmer. Currently using mmseq2 
#     # # FIrst run MUMmer to get ANI file before running mobmess (otherwise mobmess raises an error when trying to run mummer)
#     # echo "Running MUMmer to get ANI file before running mobmess..."
#     # nucmer --maxmatch --minmatch=16 -p $results_dir/mummer_align -t 40 $input_dir/predicted_plasmids_redundant.fa $input_dir/predicted_plasmids_redundant.fa

#     # echo "Running show-coords to parse delta file..."
#     # # Then parse the file
#     # show-coords -c -l -T  $results_dir/mummer_align.delta > $results_dir/mummer_align_stats.tab
#     # # This is how mummer_align_stats.tab looks like (first in TAGS is Ref, second is Query):
#     # # [S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[COV R]	[COV Q]	[TAGS]
#     # # 1	    2212	1   	2212	2212	2212	100.00	2212	2212	100.00	100.00	L1_2C_000000001932	L1_2C_000000001932
#     # # Get the necessary fields for mobmess, and then run the python sort script to join local alignments for the same plasmids together and have the same input as mobmess expects 
#     # # This is how it should look like:   Reference_Sequence	Query_Sequence	C	I_local	I_global
#     # echo "Getting pairwise alignments into ANI file..."
#     # cat $results_dir/mummer_align_stats.tab | python3 /mnt/synology/ARNAU/SCRIPTS/PLASMIDOME/sort_pairwise_alignments.py > $results_dir/predicted_plasmids_pairwise_alignments.txt

#     mkdir -p $results_dir/mmseq2_database

#     echo "Creating mmseqs database..."
#     mmseqs createdb -v 2 $input_dir/predicted_plasmids_redundant.fa $results_dir/mmseq2_database/predicted_plasmids_DB

#     echo "Creating mmseqs index..."
#     mmseqs createindex --search-type 3 -v 2 $results_dir/mmseq2_database/predicted_plasmids_DB /tmp

#     # Run the search of the plasmid list against itself
#     echo "Running mmseqs search..."
#     mmseqs search -a --alignment-mode 3 -v 2 $results_dir/mmseq2_database/predicted_plasmids_DB $results_dir/mmseq2_database/predicted_plasmids_DB $results_dir/mmseq2_database/alignment_DB /tmp

#     echo "Running mmseqs convertalis to obtain columns: query,target,fident(fraction of identity),qstart,qend,tstart,tend,qlen,tlen,cigar,evalue,qcov,tcov"
#     mmseqs convertalis --search-type 3 -v 2 --format-mode 4 --format-output "query,target,fident,qstart,qend,tstart,tend,qlen,tlen,cigar,evalue,qcov,tcov" \
#       --threads $threads $results_dir/mmseq2_database/predicted_plasmids_DB $results_dir/mmseq2_database/predicted_plasmids_DB $results_dir/mmseq2_database/alignment_DB $results_dir/alignment_stats.tsv

#     echo "Sorting pairwise alignments..."
#     cat $results_dir/alignment_stats.tsv | python3 /mnt/synology/ARNAU/SCRIPTS/PLASMIDOME/sort_pairwise_alignments_mmseq2.py > $results_dir/predicted_plasmids_pairwise_alignments.txt

    
# fi  

if [ "$step" == "2" ]; then

    echo "Running mobmess systems with precomputed MUMmer alignments..."
    mobmess systems \
    --sequences $input_dir/clean_plasmids_redundant.fa \
    --complete $results_dir/plasmids_circular.txt \
    --output $results_dir/plasmids-mobmess \
    --threads $threads \
    --ani $results_dir/clean_redundant_plasmids_pairwise_alignments.txt

fi

if [ "$step" == "3" ]; then

  echo "Getting deduplicated list of plasmids..."
  # Get the list of representative plasmids of each cluster
  # But, if the plasmid type is a fragment and it's part of a system (so there are complete plasmids there), ignore it 
  awk '$1 == $5 && (($6 != "" && $3 != "fragment") || $6 == "") {print $5}' $results_dir/plasmids-mobmess_clusters.txt > $input_dir/nonredundant_plasmids_list.txt

  # Also get the deduplicated list of only COMPLETE plasmids (not fragment)
  awk '$1 == $5 && ($3 != "fragment") {print $5}' $results_dir/plasmids-mobmess_clusters.txt > $input_dir/nonredundant_complete_plasmids_list.txt

  echo "Filtering plasmids from the redundant plasmids file..."
  grep --no-group-separator -A 1 -f $input_dir/nonredundant_plasmids_list.txt $input_dir/clean_plasmids_redundant.fa > $input_dir/nonredundant_plasmids.fa
  grep --no-group-separator -A 1 -f $input_dir/nonredundant_complete_plasmids_list.txt $input_dir/clean_plasmids_redundant.fa > $input_dir/nonredundant_complete_plasmids.fa

  echo "Filtering the plasmid info files..."
  head -n 1 $input_dir/clean_redundant_plasmids_info.txt > $input_dir/nonredundant_plasmids_info.txt
  grep -f $input_dir/nonredundant_plasmids_list.txt $input_dir/clean_redundant_plasmids_info.txt >> $input_dir/nonredundant_plasmids_info.txt


fi