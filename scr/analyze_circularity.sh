#!/bin/bash

#######WARNING######
# This sctipt neds the conda environment bowtie2 activated

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -S <step> -P <passwd> ]"
    echo
    echo "Options:"
    echo "  -i  INPUT   Input directory (with the KNEAD DATA reads for each sample)"
    echo "              Default: BASE_DIR/SYMLINKS_KNEADDATA/"
    echo "  -o  OUTPUT  Output directory to store results. Default: BASE_DIR/CIRCULARITY/"
    echo "  -p  PLASMIDS_DIRECTORY   Directory where the predicted plasmids are stored (default: BASE_DIR/PLASMIDS/)"
    echo "  -t  THREADS Number of threads (default: 1)"
    echo "  -s  SAMPLE LIST "
    echo "              File with a list of the samples to run this script for. "
    echo "              Default: BASE_DIR/SAMPLE_LISTS/POP.txt"
    echo "  -r  READS DIR   TEMP directory where to store ordered reads (default: BASE_DIR/SYMLINKS_KNEADDATA/tmp/POP)"
    echo "  -P  PSSWD   Password to run sudo commands"
    echo "  -v  ViralVerify directory     Directory with the results of ViralVerify (default: BASE_DIR/VIRALVERIFY)"
    echo "  -a  ALIGNMENT RESULTS    File with the results of the mmseq2 alignments (default: BASE_DIR/MOBMESS/alignment_stats.tsv)"
    echo "  -e  STEP    Step to run. Default: all. Options:"
    echo "                  - 1: Align reads to predicted plasmids of the same sample"
    echo "                  - 2: Filter SAM file to obtain mapping of reverse-forward reads"
    echo "                  - 3: Compute coverage of the plasmids"
    echo "                  - 4: Determine circularity of the plasmids and annotate it to info file in the plasmids directory"
    echo "  -h          Display this help message"
}

# Get options to run different parts of the strainphlan pipeline
while getopts i:o:t:p:s:P:v:a:e:r:h flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) results_dir=${OPTARG};;
        p) plasmids_dir=${OPTARG};;
        t) threads=${OPTARG};;
        s) sample_list=${OPTARG};;
        P) passwd=${OPTARG};;
        e) step=${OPTARG};;
        r) reads_dir=${OPTARG};;
        v) viralverify_dir=${OPTARG};;
        a) alignment_results=${OPTARG};;
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
  input_dir="BASE_DIR/SYMLINKS_KNEADDATA"
fi

if [ -z "$results_dir" ]; then
  results_dir="BASE_DIR/CIRCULARITY"
fi

if [ -z "$sample_list" ]; then
  sample_list="BASE_DIR/SAMPLE_LISTS/POP.txt"
fi

if [ -z "$plasmids_dir" ]; then
  plasmids_dir="BASE_DIR/PLASMIDS/"
fi

if [ -z "$step" ]; then
  step="all"
fi

if [ -z "$reads_dir" ]; then
  reads_dir="BASE_DIR/SYMLINKS_KNEADDATA/tmp/POP"
fi

if [ -z "$viralverify_dir" ]; then
  viralverify_dir="BASE_DIR/VIRALVERIFY"
fi

if [ -z "$alignment_results" ]; then
  alignment_results="BASE_DIR/MOBMESS/alignment_stats.tsv"
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
    if [ "$step" == "1" ] || [ "$step" == "all" ]; then
        echo "Analyzing sample $sample"
        if ! [ -d "$results_dir/${sample}/" ]; then
            mkdir -pv "$results_dir/${sample}"
        else
            echo "$results_dir/${sample}/ already exists, recalculating it..."
            #continue
        fi

        # Check if the .1.bt2 file exists
        if [ ! -f "$plasmids_dir/${sample}/${sample}-predicted-plasmids.1.bt2" ]; then
            echo "Indexing plasmids for sample $sample"
            bowtie2-build "$plasmids_dir/${sample}/${sample}-predicted-plasmids.fa" "$plasmids_dir/${sample}/${sample}-predicted-plasmids"
        else
            echo "Index already exists for sample $sample"
        fi

        echo "Sorting reads for sample $sample"
        # Sort reads by name so paired reads are in the same order (for bowte2 to perform well)
        zcat "$input_dir/${sample}/${sample}_kneaddata_paired_1.fastq.gz" | paste - - - - | sort -k1,1 -t " " --temporary-directory=$reads_dir | tr "\t" "\n" > "$reads_dir/${sample}_kneaddata_paired_1-sorted.fastq" &
        pid1=$!

        zcat "$input_dir/${sample}/${sample}_kneaddata_paired_2.fastq.gz" | paste - - - - | sort -k1,1 -t " " --temporary-directory=$reads_dir | tr "\t" "\n" > "$reads_dir/${sample}_kneaddata_paired_2-sorted.fastq" &
        pid2=$!

        # Wait for both background tasks to complete
        wait $pid1
        wait $pid2

        echo "Aligning sample $sample"
        # run alignment
        time bowtie2 -x "$plasmids_dir/${sample}/${sample}-predicted-plasmids" \
        -1 "$reads_dir/${sample}_kneaddata_paired_1-sorted.fastq" \
        -2 "$reads_dir/${sample}_kneaddata_paired_2-sorted.fastq" \
        -S "$results_dir/${sample}/${sample}.sam" \
        -p "$threads" \
        --no-unal > "$results_dir/${sample}/${sample}.log"  2>&1

        # Clear swap
        echo 'metalab00' | sudo -S swapoff -a; sudo -S swapon -a

        # Remove the sorted reads
        rm "$reads_dir/${sample}_kneaddata_paired_1-sorted.fastq" "$reads_dir/${sample}_kneaddata_paired_2-sorted.fastq"
    fi

    if [ "$step" == "2" ] || [ "$step" == "all" ]; then
        echo "Filtering SAM file for sample $sample"

        # Sort the SAM file 
        samtools sort --output-fmt=SAM "$results_dir/${sample}/${sample}.sam" > "$results_dir/${sample}/${sample}-sorted.sam"

        # Store the sam header
        samtools view --header-only "$results_dir/${sample}/${sample}-sorted.sam" > "$results_dir/${sample}/${sample}-header.sam"

        # Filter SAM file to obtain mapping of reverse-forward reads
        # samtools to retrieve reads with mapping quality > 20 (0.99 probability of good alignment)
        # awk to retrieve reads that: -are both mapped  - one in reverse one in forward  - incorrect insert size and also ensure both reads are in the same chromosome (same plasmid) (rnext == "=")
        # then retrieve reads that are oriented in reverse-forward orientation (POS < PNEXT if reverse read - 81 and 145; POS > PNEXT if forward read - 97 and 161) 
        samtools view -q 20  "$results_dir/${sample}/${sample}-sorted.sam" --threads "$threads" | awk '{
                                                                                                    if (($2 == 97 || $2 == 145 || $2 == 81 || $2 == 161) && ($7 == "=")) {
                                                                                                        if (($2 == 81 || $2 == 145) && ($4 < $8)) {
                                                                                                            print $0
                                                                                                        } else if (($2 == 97 || $2 == 161) && ($4 > $8)) {
                                                                                                            print $0
                                                                                                        }
                                                                                                    }
                                                                                                }' > "$results_dir/${sample}/${sample}-RF.tmp.sam"

        # Concatenate header and Reverse-Forward (RF) reads
        cat "$results_dir/${sample}/${sample}-header.sam" "$results_dir/${sample}/${sample}-RF.tmp.sam" > "$results_dir/${sample}/${sample}-RF.sam"

        # Remove the temporary file
        rm "$results_dir/${sample}/${sample}-RF.tmp.sam"
    fi

    if [ "$step" == "3" ] || [ "$step" == "all" ]; then
    
        #Calculate plasmids coverage
        echo "Calculating coverage for sample $sample"
        samtools coverage "$results_dir/${sample}/${sample}-sorted.sam" > "$results_dir/${sample}/${sample}-coverage.txt"
        
    fi

    if [ "$step" == "4" ] || [ "$step" == "all" ]; then

        echo "Analyzing circularity for sample $sample"

        # prepare a plasmid information file to store the results
        echo -e "Plasmid_ID\tPlasX_Score\tLength\tCoverage\tRF_pairs\tCircularity" > "$plasmids_dir/${sample}/${sample}-plasmid-info.txt"

        # iterate for each plasmid in the sample
        for plasmidID in $(awk 'NR > 1 {print $1}' "$results_dir/${sample}/${sample}-coverage.txt")
        do
            # get the plasX prediction value from the file
            plasX_score=$(awk '$1 == "'$plasmidID'"{print $2}' "$plasmids_dir/${sample}/${sample}-predicted-plasmids.txt")

            # get the coverage of the plasmid
            coverage=$(awk '$1 == "'$plasmidID'" {print $7}' "$results_dir/${sample}/${sample}-coverage.txt")
            
            # get the length of the plasmid
            length=$(awk '$1 == "'$plasmidID'" {print $3}' "$results_dir/${sample}/${sample}-coverage.txt")

            # calculate the median insert size of inward oriented pairs: select reads aligned to the plasmid in the sam file. 
            #Then with awk, select those where $7 (the plasmid of the pair) is the same: reads aligned to the same plasmid
            # Then use samtools stats to compute stats, whic include insert size and that can be selected with | grep ^IS | cut -f 2-
            # And then another awk to calculate median insert size: iterate per line, and store the insert sizes in an array ($3 is number of inward oriented pairs with the insert size $1)
            # Finally calculate median at the end of the awk command
            medianIS=$(grep "$plasmidID" "$results_dir/${sample}/${sample}-sorted.sam" | awk 'NR == 1 {print $0} $7 == "=" {print $0}' | samtools stats | grep ^IS | cut -f 2- | LC_NUMERIC=C awk '{
                                                                                                        for (j = 0; j < $3; j++) {
                                                                                                            a[i++] = $1;
                                                                                                        }
                                                                                                    }
                                                                                                    END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')

            # Calculate D: contig length - 3 times the median insert size
            D=$(echo "$length - 3 * $medianIS" | bc)

            #filter the Reverse-Forward reads that map to the plasmid, and use awk to check that the gap ($8 PNEXT - $4 POS) is greater than D
            RFcount=$(samtools view "$results_dir/${sample}/${sample}-RF.sam" | grep "$plasmidID" | awk 'BEGIN {count=0}{diff = $8 - $4; if(diff >= '$D') {count++;}}END {print count;}')

            # Determine circularity if at least one pair of reads is separated by a distance greater than D
            if [ "$RFcount" -gt 0 ]; then
                circularity="circular"
            else
                circularity=""
            fi

            echo -e "$plasmidID\t$plasX_score\t$length\t$coverage\t$RFcount\t$circularity" >> "$plasmids_dir/${sample}/${sample}-plasmid-info.txt"

        done

    fi


    
    

done


