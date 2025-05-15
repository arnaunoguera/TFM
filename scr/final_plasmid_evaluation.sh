#!/bin/bash


# Script to: run blastn alignments of each plasmid against the PLSDB database to check if certain contigs are plasmids for sure
#######WARNING######
# This sctipt neds the conda environment blast activated

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -S <step> -P <passwd> ]"
    echo
    echo "Options:"
    echo "  -p  PLASMIDS_DOIR   Input dir (with redundant plasmids fasta file)"
    echo "                      Default: BASE_DIR/PLASMIDS"
    echo "  -P  PLSDB DIR       Directory with results of mapping against PLSDB. Default: BASE_DIR/PLSDB_BLAST"
    echo "  -v  ViralVerify DIR Directory with results of ViralVerify. Default: BASE_DIR/VIRALVERIFY"
    echo "  -b  BACTERIAL CHR DIR   Directory with results of mapping against bacterial chromosomes. Default: BASE_DIR/BACTERIAL_CHROMOSOME_FILTERING"
    echo "  -m  MOB-SUITE DIR   Directory with results of mob-suite. Default: BASE_DIR/MOB_SUITE"
    echo "  -s  STEP            Step to run. Default: 1. Options:"
    echo "                      - 1: Split the plasmid info file into n parts to parallelize"  
    echo "                      - 2: Unite all the information from all sources into a plasmid info file"
    echo "                      - 3: Merge all the information into a single file"
    echo "                      - 4: Get a list of plasmid-related cog and pfam functions to add"
    echo "  -t  PARTITION    Partition number for multithreading (Default: 00)"
    echo "  -n  N_PARTITIONS   Number of partitions to split the input file (default: 30)"
    echo "  -h          Display this help message"
}

# Get options to run different parts of the strainphlan pipeline
while getopts p:P:v:b:m:s:t:n:h flag
do
    case "${flag}" in
        p) plasmid_dir=${OPTARG};;
        P) PLSDB_dir=${OPTARG};;
        v) viral_dir=${OPTARG};;
        b) bacterial_dir=${OPTARG};;
        m) MOB_dir=${OPTARG};;
        s) step=${OPTARG};;
        t) partition=${OPTARG};;
        n) n_partitions=${OPTARG};;
        h) show_help
           exit 0;;
        *) show_help
           exit 1;;
    esac
done


###### DEFINE VARIABLES TO RUN PLASX ######


if [ -z "$plasmid_dir" ]; then
  plasmid_dir="BASE_DIR/PLASMIDS"
fi

if [ -z "$PLSDB_dir" ]; then
  PLSDB_dir="BASE_DIR/PLSDB_BLAST"
fi

if [ -z "$viral_dir" ]; then
  viral_dir="BASE_DIR/VIRALVERIFY"
fi

if [ -z "$bacterial_dir" ]; then
  bacterial_dir="BASE_DIR/BACTERIAL_CHROMOSOME_FILTERING"
fi

if [ -z "$MOB_dir" ]; then
  MOB_dir="BASE_DIR/MOB_SUITE"
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

#### MAIN SCRIPT 

if [ "$step" == "1" ]; then

    echo "Splitting input file into $n_partitions parts to parallelize..."
    n_plasmids=$(wc -l $plasmid_dir/redundant_plasmid_joint_info.txt | cut -d " " -f 1)
    n_per_partition=$(( n_plasmids / n_partitions))
    split -da 2 -l $n_per_partition $plasmid_dir/redundant_plasmid_joint_info.txt $plasmid_dir/parallel_runs/redundant_plasmid_joint_info_ --additional-suffix=".txt" --verbose 

fi

if [ "$step" == "2" ]; then
    
    n=0
    # Iterate on each line of the plasmid info file
    while read -r line; do

        if [ $partition == "00" ] && [ $n -eq 0 ]; then # print the header directly
            n=1
            # The header is Plasmid_ID\tPlasX_Score\tLength\tCoverage\tRF_pairs\tCircularity\tCircularity_Source
            # Now add: ViralVerify_prediction\tViralVerify_Score\tBacterial_Chr\tBacterial_Chr_cov\tBacterial_Chr_ident\tPLSDB_hit\tPLSDB_cov\tPLSDB_ident
            echo -e "$line\tViralVerify_prediction\tViralVerify_Score\tBacterial_Chr\tBacterial_Chr_cov\tBacterial_Chr_ident\tPLSDB_hit\tPLSDB_cov\tPLSDB_ident\trep_type(s)\trep_type_accession(s)\trelaxase_type(s)\trelaxase_type_accession(s)\tmpf_type\tmpf_type_accession(s)\torit_type(s)\torit_accession(s)\tpredicted_mobility" > $plasmid_dir/parallel_runs/redundant_plasmid_final_info_$partition.txt
        else
            if [ $((n % 100)) -eq 0 ]; then
                echo "Processing plasmid $n"
            fi
            n=$((n+1)) # update counter

            plasmidID=$(echo -e "$line" | cut -f 1)
            
            # get prediction column and viralverify score column from the viralverify results
            viralverify_info=$(grep "$plasmidID" $viral_dir/predicted_plasmids_redundant_result_table.csv | cut -f 2,5 -d ",")
            viralverify_prediction=$(echo -e "$viralverify_info" | cut -f 1 -d ",")
            viralverify_score=$(echo -e "$viralverify_info" | cut -f 2 -d ",")

            #get the topmost hit (in coverage) from the bacterial chromosome for the plasmid
            bacterial_decontamination=$(grep "$plasmidID" $bacterial_dir/blast_results_COV90.tsv | sort -k 4 -g -r | head -n 1)

            # check if there were any results 
            if [ -z "$bacterial_decontamination" ]; then
                bacterial_cov="NA"
                bacterial_ident="NA"
                bacterial_hits="NA"
            else 
                bacterial_cov=$(echo -e "$bacterial_decontamination" | cut -f 4)
                bacterial_ident=$(echo -e "$bacterial_decontamination" | cut -f 3)
                bacterial_hits=$(grep "$plasmidID" $bacterial_dir/bacterial_chromosome_hits.tsv | cut -f 2 | paste -sd ',')
            fi

            # Next, check the plasmid results for the PLSDB mapping 
            PLSDB_info=$(grep "$plasmidID" $PLSDB_dir/blast_results_COV90.tsv | sort -k 4 -g -r | head -n 1)
            # check if there were any results
            if [ -z "$PLSDB_info" ]; then
                PLSDB_cov="NA"
                PLSDB_ident="NA"
                PLSDB_hits="NA"
            else
                PLSDB_cov=$(echo -e "$PLSDB_info" | cut -f 4)
                PLSDB_ident=$(echo -e "$PLSDB_info" | cut -f 3)
                PLSDB_hits=$(grep "$plasmidID" $PLSDB_dir/PLSDB_hits.tsv | cut -f 2 | paste -sd ',')
            fi

            # Finally, get the mob-typer information. I want the columns 6-14:
            # rep_type(s)	rep_type_accession(s)	relaxase_type(s)	relaxase_type_accession(s)
            # mpf_type	mpf_type_accession(s)	orit_type(s)	orit_accession(s)	predicted_mobility
            # and replace - by NAs
            mobtyper_info=$(grep "$plasmidID" $MOB_dir/mobtyper_results.txt | cut -f 6-14 | sed 's/-\t/NA\t/g' )

            # Now add the information to the plasmid info file
            echo -e "$line\t$viralverify_prediction\t$viralverify_score\t$bacterial_hits\t$bacterial_cov\t$bacterial_ident\t$PLSDB_hits\t$PLSDB_cov\t$PLSDB_ident\t$mobtyper_info" >> $plasmid_dir/parallel_runs/redundant_plasmid_final_info_$partition.txt

        fi

    done < "$plasmid_dir/parallel_runs/redundant_plasmid_joint_info_$partition.txt"
fi

if [ "$step" == "3" ]; then

    echo "Merging parallelized results..."
    cat $plasmid_dir/parallel_runs/redundant_plasmid_final_info_*.txt > $plasmid_dir/redundant_plasmid_final_info.txt
    echo "Results merged in $plasmid_dir/redundant_plasmid_final_info.txt, removing temp files"
    rm $plasmid_dir/parallel_runs/redundant_plasmid_final_info_*.txt

fi

if [ "$step" == "4" ]; then

  echo "Looking for plasmid-related COG and pfam functions"
  head -n 1 $plasmid_dir/redundant_plasmids_cogs-and-pfams.txt > $plasmid_dir/plasmid_related_cogs_and_pfams.txt
  grep -i  "plasmid"  $plasmid_dir/redundant_plasmids_cogs-and-pfams.txt | grep -iv "phage" >> $plasmid_dir/plasmid_related_cogs_and_pfams.txt

  cut -f 1 $plasmid_dir/plasmid_related_cogs_and_pfams.txt > $plasmid_dir/plasmid_related_cogs_and_pfams_list.txt

  echo "Looking for plasmid-related gene calls"
  awk 'FNR==NR { ids[$1]; next } $1 in ids {print $0}' $plasmid_dir/plasmid_related_cogs_and_pfams_list.txt $plasmid_dir/redundant_plasmids_gene-calls.txt > $plasmid_dir/plasmid_related_gene_calls.txt

  rm $plasmid_dir/plasmid_related_cogs_and_pfams_list.txt

fi