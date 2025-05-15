#!/bin/bash


# Script to: make the pariwise alignments necessary for mobmess
# Join the circularity annotations of 3 sources: RF read pairs, ViralVerify and pairwise alignments
#######WARNING######
# This sctipt neds the conda environment mmseq2 activated

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -S <step> -P <passwd> ]"
    echo
    echo "Options:"
    echo "  -i  INPUT   Input dir (with redundant plasmids fasta file)"
    echo "              Default: BASE_DIR/PLASMIDS"
    echo "  -o  OUTPUT  Output directory to store the pairwise alignment results. Default: BASE_DIR/MOBMESS"
    echo "  -t  THREADS Number of threads (default: 1)"
    echo "  -v  ViralVerify directory     Directory with the results of ViralVerify (default: BASE_DIR/VIRALVERIFY)"
    echo "  -s  STEP    Step to run. Default: 1. Options:"
    echo "                  - 1: Run the pairwise alignments with mmseq2"
    echo "                  - 2: Split input file into 20 parts to parallelize"
    echo "                  - 3: Annotate circularit based on ViralVerify and pairwise alignments results"
    echo "                  - 4: Merge parallelized results"
    echo "  -m  PARTITION    Partition number for multithreading"
    echo "  -n  N_PARTITIONS   Number of partitions to split the input file (default: 30)"
    echo "  -h          Display this help message"
}

# Get options to run different parts of the strainphlan pipeline
while getopts i:o:t:v:n:a:s:m:r:h flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) results_dir=${OPTARG};;
        t) threads=${OPTARG};;
        s) step=${OPTARG};;
        v) viralverify_dir=${OPTARG};;
        m) partition=${OPTARG};;
        n) n_partitions=${OPTARG};;
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
  results_dir="BASE_DIR/MOBMESS"
fi

if [ -z "$step" ]; then
  step="1"
fi

if [ -z "$partition" ]; then
  partition="00"
fi

if [ -z "$n_partitions" ]; then
  n_partitions=30
fi


if [ -z "$viralverify_dir" ]; then
  viralverify_dir="BASE_DIR/VIRALVERIFY"
fi


#####MAIN SCRIPT######

if [ "$step" == "1" ]; then

    # This commented chunk was to run alignment with nucmer. Currently using mmseq2 
    # # FIrst run MUMmer to get ANI file before running mobmess (otherwise mobmess raises an error when trying to run mummer)
    # echo "Running MUMmer to get ANI file before running mobmess..."
    # nucmer --maxmatch --minmatch=16 -p $results_dir/mummer_align -t 40 $input_dir/predicted_plasmids_redundant.fa $input_dir/predicted_plasmids_redundant.fa

    # echo "Running show-coords to parse delta file..."
    # # Then parse the file
    # show-coords -c -l -T  $results_dir/mummer_align.delta > $results_dir/mummer_align_stats.tab
    # # This is how mummer_align_stats.tab looks like (first in TAGS is Ref, second is Query):
    # # [S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[COV R]	[COV Q]	[TAGS]
    # # 1	    2212	1   	2212	2212	2212	100.00	2212	2212	100.00	100.00	L1_2C_000000001932	L1_2C_000000001932
    # # Get the necessary fields for mobmess, and then run the python sort script to join local alignments for the same plasmids together and have the same input as mobmess expects 
    # # This is how it should look like:   Reference_Sequence	Query_Sequence	C	I_local	I_global
    # echo "Getting pairwise alignments into ANI file..."
    # cat $results_dir/mummer_align_stats.tab | python3 /mnt/synology/ARNAU/SCRIPTS/PLASMIDOME/sort_pairwise_alignments.py > $results_dir/predicted_plasmids_pairwise_alignments.txt

    mkdir -p $results_dir/mmseq2_database

    echo "Creating mmseqs database..."
    mmseqs createdb -v 2 $input_dir/predicted_plasmids_redundant.fa $results_dir/mmseq2_database/predicted_plasmids_DB

    echo "Creating mmseqs index..."
    mmseqs createindex --search-type 3 -v 2 $results_dir/mmseq2_database/predicted_plasmids_DB /tmp

    # Run the search of the plasmid list against itself
    echo "Running mmseqs search..."
    mmseqs search -a --alignment-mode 3 -v 2 $results_dir/mmseq2_database/predicted_plasmids_DB $results_dir/mmseq2_database/predicted_plasmids_DB $results_dir/mmseq2_database/alignment_DB /tmp

    echo "Running mmseqs convertalis to obtain columns: query,target,fident(fraction of identity),qstart,qend,tstart,tend,qlen,tlen,cigar,evalue,qcov,tcov"
    mmseqs convertalis --search-type 3 -v 2 --format-mode 4 --format-output "query,target,fident,qstart,qend,tstart,tend,qlen,tlen,cigar,evalue,qcov,tcov" \
      --threads $threads $results_dir/mmseq2_database/predicted_plasmids_DB $results_dir/mmseq2_database/predicted_plasmids_DB $results_dir/mmseq2_database/alignment_DB $results_dir/alignment_stats.tsv

    echo "Sorting pairwise alignments..."
    cat $results_dir/alignment_stats.tsv | python3 /mnt/synology/ARNAU/SCRIPTS/PLASMIDOME/sort_pairwise_alignments_mmseq2.py > $results_dir/predicted_plasmids_pairwise_alignments.txt

fi

if [ "$step" == "2" ]; then
  echo "Splitting input file into $n_partitions parts to parallelize..."
  mkdir -p $input_dir/parallel_runs
  n_plasmids=$(wc -l $input_dir/redundant_plasmid_info.txt | cut -d " " -f 1)
  n_per_partition=$(( (n_plasmids - 1 ) / n_partitions))
  split -da 2 -l $n_per_partition $input_dir/redundant_plasmid_info.txt $input_dir/parallel_runs/redundant_plasmid_info_ --additional-suffix=".txt" --verbose 
fi


if [ "$step" == "3" ]; then
  echo "Annotating circularity based on ViralVerify and pairwise alignments results..."
  n=0
  # Create a file like the current plasmid info file but with all circulairty sources
  # Read lines of the plasmid info file 
  while read -r line; do
    if [ $partition == "00" ] && [ $n -eq 0 ]; then # print the header directly
      n=1
      echo -e "$line\tCircularity_Source" > $input_dir/parallel_runs/redundant_plasmid_joint_info_$partition.txt
    else
      if [ $((n % 100)) -eq 0 ]; then
        echo "Processing plasmid $n"
      fi
      n=$((n+1)) # update counter
 
      # The header is Plasmid_ID\tPlasX_Score\tLength\tCoverage\tRF_pairs\tCircularity\tCircularity_Source
      # Extract circularity and if it was already annotated as circular, just print the line
      circularity=$(echo -e "$line" | cut -f 6)
      if [ "$circularity" == "circular" ]; then
        echo -e "$line\tRF pairs" >> $input_dir/parallel_runs/redundant_plasmid_joint_info_$partition.txt
      else
        plasmidID=$(echo -e "$line" | cut -f 1)
        PlasX_score=$(echo -e "$line" | cut -f 2)
        Length=$(echo -e "$line" | cut -f 3)
        Coverage=$(echo -e "$line" | cut -f 4)
        RF_pairs=$(echo -e "$line" | cut -f 5)

        #Check if the ends of the plasmid are overlapping
        circ_self=$(grep "^$plasmidID" $results_dir/alignment_stats.tsv | awk -v id="$plasmidID" '
            BEGIN {counter = 0}
            $1 == id && $2 == id && (($4 == 1 && $5 != $8 && $6 != 1 && $7 == $8) || ($4 != 1 && $5 == $8 && $6 == 1 && $7 != $8)) {
                counter ++
            }
            END {print counter}
        ' )
        if [ $circ_self -gt 0 ]; then
            circularity="circular"
            circularity_source="overlapping ends"
            echo "Circular plasmid (Self alignment): $plasmidID" 
        else 
          # We know the plasmid isn't circular according to RF pairs
          # check out the alignments of the plasmid as query ($1) in mmseq2 alignments, to see if the start and end coordinates align contiguously in another plasmid
          # file header: query	target	fident	qstart	qend	tstart	tend	qlen	tlen	cigar	evalue	qcov	tcov
          # get alignments that start or end in the first or last position of the query 
          #awk -v id="$plasmidID" '$1 == id && $2 != id && ($4 == 1 || $4 == $8 || $5 == 1 || $5 == $8) {print $0}'

          # Use python script detect_alignment_circularity.py to determine if the plasmid is circular based on the alignments
          circ_other_plasmids=$(grep "^$plasmidID" $results_dir/alignment_stats.tsv | python /mnt/synology/ARNAU/SCRIPTS/PLASMIDOME/detect_alignment_circularity.py)

          # circ_other_plasmids=$(grep "^$plasmidID" $results_dir/alignment_stats.tsv | awk -v id="$plasmidID" '
          #     BEGIN {counter = 0}
          #     $1 == id && $2 != id && ($4 == 1 || $4 == $8 || $5 == 1 || $5 == $8) {
          #         if (($2 SUBSEP "1") in prev) {
          #             if (($6 == prev[$2 SUBSEP "1"] + 1 || $6 == prev[$2 SUBSEP "1"] - 1) && ($7 == prev[$2 SUBSEP "2"] + 1 || $7 == prev[$2 SUBSEP "2"] - 1)) {
          #                 counter++
          #             }
          #         }
          #         prev[$2 SUBSEP "1"] = $6
          #         prev[$2 SUBSEP "2"] = $7
          #     }
          #     END {print counter}
          # ' )
          if [ $circ_other_plasmids == "1" ]; then
            circularity="circular"
            circularity_source="alignment to other plasmids"
            echo "Circular plasmid (Pairwise alignment): $plasmidID"
          else
            # Also determine circularity if the plasmid is detected as circular by ViralVerify (+ in column 4)
            viralverify_circularity=$(grep "$plasmidID" $viralverify_dir/predicted_plasmids_redundant_result_table.csv | cut -d "," -f 4)
            if [ "$viralverify_circularity" == "+" ]; then
              circularity="circular"
              circularity_source="ViralVerify"
              echo "Circular plasmid (ViralVerify): $plasmidID"
            else
              circularity=""
              circularity_source="NA"
            fi
          fi
        fi
        echo -e "$plasmidID\t$PlasX_score\t$Length\t$Coverage\t$RF_pairs\t$circularity\t$circularity_source" >> $input_dir/parallel_runs/redundant_plasmid_joint_info_$partition.txt
      fi
    fi
  done < "$input_dir/parallel_runs/redundant_plasmid_info_$partition.txt"

  echo "Results for partition $partition are ready in file $input_dir/parallel_runs/redundant_plasmid_joint_info_$partition.txt"

fi  

if [ "$step" == "4" ]; then
  echo "Merging parallelized results..."
  cat $input_dir/parallel_runs/redundant_plasmid_joint_info_*.txt > $input_dir/redundant_plasmid_joint_info.txt
  echo "Results merged in $input_dir/redundant_plasmid_joint_info.txt, removing temp files"
  rm $input_dir/parallel_runs/redundant_plasmid_joint_info_*.txt
fi