#!/bin/bash

#######WARNING######
# This sctipt neds the conda environment bowtie2 activated

# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -A <bam_dir> -c ]"
    echo
    echo "Options:"
    echo "  -i  INPUT   Input dir (file with KNEAD DATA reads, with a folder per sample)"
    echo "              Default: BASE_DIR/SYMLINKS_KNEADDATA"
    echo "  -p  PLASMIDS DIR Directory with a fasta file with all the reference (non-redunddant) plasmids to align the reads to"
    echo "              Default: BASE_DIR/PLASMIDS"
    echo "  -o  OUTPUT  Output directory to store results. Default: BASE_DIR/PROFILING"
    echo "  -A  ALTERNATIVE OUTPUT Alternative output directory to store BAM files. Default: same as -o"
    echo "  -t  THREADS Number of threads (default: 1)"
    echo "  -w  WHICH PLASMIDS  For which plasmids to run. Options: all, complete, both. Default: both"
    echo "  -s  SAMPLE LIST "
    echo "              File with a list of the samples to run this script for. "
    echo "              Default: BASE_DIR/SAMPLE_LISTS/POP.txt"
    echo "  -e  STEP    Step to run. Default: 1. Options:"
    echo "                  - 1: Align reads to all prediced plasmids and compute the coverage"
    echo "                  - 2: Merge the coverage files of all samples (optional -- not run the final time )"
    echo "                  - 3: Merge the coverage but for ALL samples of all 3 folders (needs samplelist file with sample ID and group name)"
    echo "                  - 4: Create a BED file with the gene annotations for gene profiling"
    echo "                  - 5: index the bam files"
    echo "                  - 6: run bedtools multicov to exytract the number of reads mapped to each gene region (requires conda env bedtools)"
    echo "  -F FINAL ANNOTATION   Directory with final plasmid annotations (gene calls). Default: BASE_DIR/ANVIO/FINAL_PLASMIDS/"
    echo "  -h          Display this help message"
}


# Get options to run different parts of the strainphlan pipeline
while getopts i:o:p:t:s:e:A:P:w:F:h flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) results_dir=${OPTARG};;
        p) plasmids_dir=${OPTARG};;
        t) threads=${OPTARG};;
        A) bam_dir=${OPTARG};;
        s) sample_list=${OPTARG};;
        e) step=${OPTARG};;
        F) final_annotation=${OPTARG};;
        w) which_plasmids=${OPTARG};;
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
  input_dir="BASE_DIR/SYMLINKS_KNEADDATA"
fi

if [ -z "$results_dir" ]; then
  results_dir="BASE_DIR/PROFILING"
fi

if [ -z "$plasmids_dir" ]; then
  plasmids_dir="BASE_DIR/PLASMIDS"
fi

if [ -z "$step" ]; then
  step="1"
fi

if [ -z "$sample_list" ]; then
  sample_list="BASE_DIR/SAMPLE_LISTS/POP.txt"
fi

if [ -z "$bam_dir" ]; then
  bam_dir="$results_dir"
fi

if [ -z "$step" ]; then
  step="1"
fi

if [ -z "$final_annotation" ]; then
  final_annotation="BASE_DIR/ANVIO/FINAL_PLASMIDS"
fi

if [ -z "$which_plasmids" ]; then
  which_plasmids="both"
fi

# Check if which plasmids is both
if [ "$which_plasmids" == "both" ]; then
  run_all=1
  run_complete=1
fi
if [ "$which_plasmids" == "all" ]; then
  run_all=1
  run_complete=0
fi
if [ "$which_plasmids" == "complete" ]; then
  run_all=0
  run_complete=1
fi

#####MAIN SCRIPT######

if [ -d $input_dir ]; then
  echo "Input directory is $input_dir"
else
  echo "ERROR: Couldn't find input directory. Exitting..."
  exit 1
fi


# index the plasmids if it's not done already, both the nonredundant list and the one with only complete plasmids 
if ! [ -d "$plasmids_dir/bowtie2_index" ]; then
    mkdir -pv "$plasmids_dir/bowtie2_index"
    bowtie2-build "$plasmids_dir/nonredundant_plasmids.fa" "$plasmids_dir/bowtie2_index/nonredundant_plasmids" --threads "$threads"
    bowtie2-build "$plasmids_dir/nonredundant_complete_plasmids.fa" "$plasmids_dir/bowtie2_index/nonredundant_complete_plasmids" --threads "$threads"
else 
    echo "Plasmids already indexed, proceeding to profiling..."
fi



if [ "$step" == "1" ] || [ "$step" == "all" ]; then

  for sample in $(cat $sample_list);
  do
    echo "Analyzing sample $sample"

    if ! [ -d "$results_dir/${sample}/" ]; then
        mkdir -pv "$results_dir/${sample}"
        if [ "$bam_dir" != "$results_dir" ]; then
            mkdir -pv "$bam_dir/${sample}"
        fi
    else
        echo "$results_dir/${sample}/ already exists, skipping to the next sample..."
        #continue
    fi

    echo "Sorting reads for sample $sample"
    # Sort reads by name so paired reads are in the same order (for bowte2 to perform well)
    zcat "$input_dir/${sample}/${sample}_kneaddata_paired_1.fastq.gz" | paste - - - - | parsort --parallel=$threads -k1,1 -t " " --temporary-directory=/tmp/ | tr "\t" "\n" > "$input_dir/${sample}_kneaddata_paired_1-sorted.fastq" 

    zcat "$input_dir/${sample}/${sample}_kneaddata_paired_2.fastq.gz" | paste - - - - | parsort --parallel=$threads -k1,1 -t " " --temporary-directory=/tmp/ | tr "\t" "\n" > "$input_dir/${sample}_kneaddata_paired_2-sorted.fastq"

    if [ $run_all -eq 1 ]; then 

      echo "Aligning sample $sample"
      # run alignment
      time bowtie2 -x "$plasmids_dir/bowtie2_index/nonredundant_plasmids" \
      -1 "$input_dir/${sample}_kneaddata_paired_1-sorted.fastq" \
      -2 "$input_dir/${sample}_kneaddata_paired_2-sorted.fastq" \
      -S "$bam_dir/${sample}/${sample}-profiling.sam" \
      -p "$threads" \
      --no-unal > "$results_dir/${sample}/${sample}.log" 2>&1

      # Sort the SAM file 
      samtools sort --output-fmt=BAM "$bam_dir/${sample}/${sample}-profiling.sam" > "$bam_dir/${sample}/${sample}-profiling-sorted.bam"

      # if results dir and sam dir are different, make a link to the sam file in the results dir
      if [ "$bam_dir" != "$results_dir" ]; then
          ln -s "$bam_dir/${sample}/${sample}-profiling.sam" "$results_dir/${sample}/${sample}-profiling.sam"
          ln -s "$bam_dir/${sample}/${sample}-profiling-sorted.bam" "$results_dir/${sample}/${sample}-profiling-sorted.bam"
      fi

      # Get the coverage of each plasmid
      echo "Getting coverage for sample $sample"
      samtools coverage -o "$results_dir/${sample}/${sample}-profiling-coverage.txt" "$results_dir/${sample}/${sample}-profiling-sorted.bam" 

    fi
    
    if [ $run_complete -eq 1 ]; then 
      # run alignment only for complete plasmids
      echo "Aligning sample $sample only complete plasmids"
      # run alignment
      time bowtie2 -x "$plasmids_dir/bowtie2_index/nonredundant_complete_plasmids" \
      -1 "$input_dir/${sample}_kneaddata_paired_1-sorted.fastq" \
      -2 "$input_dir/${sample}_kneaddata_paired_2-sorted.fastq" \
      -S "$bam_dir/${sample}/${sample}-profiling-only-complete.sam" \
      -p "$threads" \
      --no-unal > "$results_dir/${sample}/${sample}-only-complete.log" 2>&1

      # Sort the SAM file 
      samtools sort --output-fmt=BAM "$bam_dir/${sample}/${sample}-profiling-only-complete.sam" > "$bam_dir/${sample}/${sample}-profiling-sorted-only-complete.bam"

      # if results dir and sam dir are different, make a link to the sam file in the results dir
      if [ "$bam_dir" != "$results_dir" ]; then
          ln -s "$bam_dir/${sample}/${sample}-profiling-only-complete.sam" "$results_dir/${sample}/${sample}-profiling-only-complete.sam"
          ln -s "$bam_dir/${sample}/${sample}-profiling-sorted-only-complete.bam" "$results_dir/${sample}/${sample}-profiling-sorted-only-complete.bam"
      fi

      # Get the coverage of each plasmid
      echo "Getting coverage for sample $sample"
      samtools coverage -o "$results_dir/${sample}/${sample}-profiling-coverage-only-complete.txt" "$results_dir/${sample}/${sample}-profiling-sorted-only-complete.bam" 

    fi

    # Remove the sorted reads
    rm "$input_dir/${sample}_kneaddata_paired_1-sorted.fastq" "$input_dir/${sample}_kneaddata_paired_2-sorted.fastq"

    # Clear swap
    echo 'metalab00' | sudo -S swapoff -a; sudo -S swapon -a 
  done
fi

if [ "$step" == "2" ] || [ "$step" == "all" ]; then

  # Define the output file
  numreads_file="$results_dir/numreads_file.tsv"
  coverage_file="$results_dir/coverage_file.tsv"
  meandepth_file="$results_dir/meandepth_file.tsv"

  # Initialize the header
  header="pCode"

  # Initialize a temporary file to store the merged columns
  temp_file_numreads=$(mktemp)
  temp_file_coverage=$(mktemp)
  temp_file_meandepth=$(mktemp)

  # Before the loop, define variable first to determine if it's the first file - to add the plasmid ID column
  first=1

  for sample in $(cat $sample_list);
  do
    echo "Merging sample $sample"

    file="$results_dir/${sample}/${sample}-profiling-coverage.txt"

    #store a sorted version of the file temporarily (remove last line, which is the header that gets put down there)
    sort $file | head -n -1 > "${file}_sorted.tmp"

    # Add the sample name to the header
    header="$header\t$sample"

    # Check if it's the first file, if so, add the plasmid ID column apart from the one of interest
    if [ $first -eq 1 ]; then
        # Extract the first column (plasmid ID) and add it to the temporary file along with the one of interest
        awk '{print $1 "\t" $4}' "${file}_sorted.tmp" > "${file}_col4.tmp"
        awk '{print $1 "\t" $6}' "${file}_sorted.tmp" > "${file}_col6.tmp"
        awk '{print $1 "\t" $7}' "${file}_sorted.tmp" > "${file}_col7.tmp"
        first=0

        # save this first file as the temporary file for each column
        mv "${file}_col4.tmp" "$temp_file_numreads"
        mv "${file}_col6.tmp" "$temp_file_coverage"
        mv "${file}_col7.tmp" "$temp_file_meandepth"

        rm "${file}_sorted.tmp" 
      
    # if it's not the first file
    else
        # Extract the 4th (numreads) column and add it to the temporary file
        awk '{print $4}' "${file}_sorted.tmp" > "${file}_col4.tmp"

        # Extract the 6th (coverage) column and add it to the temporary file
        awk '{print $6}' "${file}_sorted.tmp" > "${file}_col6.tmp"

        # Extract the 7th (meandepth) column and add it to the temporary file
        awk '{print $7}' "${file}_sorted.tmp" > "${file}_col7.tmp"

        # Paste the columns together
        paste "$temp_file_numreads" "${file}_col4.tmp" > "${temp_file_numreads}_new"
        mv "${temp_file_numreads}_new" "$temp_file_numreads"

        paste "$temp_file_coverage" "${file}_col6.tmp" > "${temp_file_coverage}_new"
        mv "${temp_file_coverage}_new" "$temp_file_coverage"

        paste "$temp_file_meandepth" "${file}_col7.tmp" > "${temp_file_meandepth}_new"
        mv "${temp_file_meandepth}_new" "$temp_file_meandepth"

        # Remove the temporary files
        rm "${file}_sorted.tmp" "${file}_col4.tmp" "${file}_col6.tmp" "${file}_col7.tmp"

    fi


  done

  # Add the header to the files
  echo -e $header > $numreads_file
  echo -e $header > $coverage_file
  echo -e $header > $meandepth_file

  # Add the content of the temporary files to the final files
  cat "$temp_file_numreads" >> $numreads_file
  cat "$temp_file_coverage" >> $coverage_file
  cat "$temp_file_meandepth" >> $meandepth_file

  # Remove the temporary files
  rm "${temp_file_numreads}" "${temp_file_coverage}" "${temp_file_meandepth}" 

fi

if [ "$step" == "3" ] || [ "$step" == "all" ]; then

  # Define the output file
  numreads_file="$results_dir/numreads.tsv"
  coverage_file="$results_dir/coverage.tsv"
  meandepth_file="$results_dir/meandepth.tsv"
  numreads_file2="$results_dir/numreads_only-complete.tsv"
  coverage_file2="$results_dir/coverage_only-complete.tsv"
  meandepth_file2="$results_dir/meandepth_only-complete.tsv"

  # Initialize the header
  header="pCode"

  # Initialize a temporary file to store the merged columns
  temp_file_numreads=$(mktemp)
  temp_file_coverage=$(mktemp)
  temp_file_meandepth=$(mktemp)
  temp_file_numreads2=$(mktemp)
  temp_file_coverage2=$(mktemp)
  temp_file_meandepth2=$(mktemp)

  # Before the loop, define variable first to determine if it's the first file - to add the plasmid ID column
  first=1

  # Iterate over each line of the file
  while IFS=$'\t' read -r sample sample_group; do
    
    # Process each line
    echo "Merging sample $sample from $sample_group"

    file="$results_dir/$sample_group/${sample}/${sample}-profiling-coverage.txt"
    file2="$results_dir/$sample_group/${sample}/${sample}-profiling-coverage-only-complete.txt"

    #store a sorted version of the file temporarily (remove first line, which is the header)
    tail -n +2 $file | sort > "${file}_sorted.tmp" &
    pid1=$!
    tail -n +2 $file2 | sort > "${file2}_sorted.tmp" &
    pid2=$!

    # Wait for both background tasks to complete
    wait $pid1
    wait $pid2

    # Add the sample name to the header
    header="$header\t$sample"

    # Check if it's the first file, if so, add the plasmid ID column apart from the one of interest
    if [ $first -eq 1 ]; then
        # Extract the first column (plasmid ID) and add it to the temporary file along with the one of interest
        awk '{print $1 "\t" $4}' "${file}_sorted.tmp" > "${file}_col4.tmp" &
        pid1=$!
        awk '{print $1 "\t" $6}' "${file}_sorted.tmp" > "${file}_col6.tmp" &
        pid2=$!
        awk '{print $1 "\t" $7}' "${file}_sorted.tmp" > "${file}_col7.tmp" &
        pid3=$!
        awk '{print $1 "\t" $4}' "${file2}_sorted.tmp" > "${file2}_col4.tmp" &
        pid4=$!
        awk '{print $1 "\t" $6}' "${file2}_sorted.tmp" > "${file2}_col6.tmp" &
        pid5=$!
        awk '{print $1 "\t" $7}' "${file2}_sorted.tmp" > "${file2}_col7.tmp" &
        pid6=$!

        first=0

        wait $pid1 $pid2 $pid3 $pid4 $pid5 $pid6

        # save this first file as the temporary file for each column
        mv "${file}_col4.tmp" "$temp_file_numreads"
        mv "${file}_col6.tmp" "$temp_file_coverage"
        mv "${file}_col7.tmp" "$temp_file_meandepth"
        mv "${file2}_col4.tmp" "$temp_file_numreads2"
        mv "${file2}_col6.tmp" "$temp_file_coverage2"
        mv "${file2}_col7.tmp" "$temp_file_meandepth2"

        rm "${file}_sorted.tmp" "${file2}_sorted.tmp" 
      
    # if it's not the first file
    else
        # Extract the 4th (numreads) column and add it to the temporary file
        awk '{print $4}' "${file}_sorted.tmp" > "${file}_col4.tmp" &
        pid1=$!
        awk '{print $4}' "${file2}_sorted.tmp" > "${file2}_col4.tmp" &
        pid2=$!

        # Extract the 6th (coverage) column and add it to the temporary file
        awk '{print $6}' "${file}_sorted.tmp" > "${file}_col6.tmp" & 
        pid3=$!
        awk '{print $6}' "${file2}_sorted.tmp" > "${file2}_col6.tmp" &
        pid4=$!

        # Extract the 7th (meandepth) column and add it to the temporary file
        awk '{print $7}' "${file}_sorted.tmp" > "${file}_col7.tmp" &
        pid5=$!
        awk '{print $7}' "${file2}_sorted.tmp" > "${file2}_col7.tmp" &
        pid6=$!

        wait $pid1 $pid2 $pid3 $pid4 $pid5 $pid6

        # Paste the columns together
        paste "$temp_file_numreads" "${file}_col4.tmp" > "${temp_file_numreads}_new" &
        pid1=$!
        paste "$temp_file_coverage" "${file}_col6.tmp" > "${temp_file_coverage}_new" &
        pid2=$!
        paste "$temp_file_meandepth" "${file}_col7.tmp" > "${temp_file_meandepth}_new" &
        pid3=$!
        paste "$temp_file_numreads2" "${file2}_col4.tmp" > "${temp_file_numreads2}_new" &
        pid4=$!
        paste "$temp_file_coverage2" "${file2}_col6.tmp" > "${temp_file_coverage2}_new" &
        pid5=$!
        paste "$temp_file_meandepth2" "${file2}_col7.tmp" > "${temp_file_meandepth2}_new" &
        pid6=$!

        wait $pid1 $pid2 $pid3 $pid4 $pid5 $pid6

        mv "${temp_file_numreads}_new" "$temp_file_numreads" &
        pid1=$!
        mv "${temp_file_coverage}_new" "$temp_file_coverage" &
        pid2=$!
        mv "${temp_file_meandepth}_new" "$temp_file_meandepth" &
        pid3=$!
        mv "${temp_file_numreads2}_new" "$temp_file_numreads2" &
        pid4=$!
        mv "${temp_file_coverage2}_new" "$temp_file_coverage2" &
        pid5=$!
        mv "${temp_file_meandepth2}_new" "$temp_file_meandepth2" &
        pid6=$!

        wait $pid1 $pid2 $pid3 $pid4 $pid5 $pid6

        # Remove the temporary files
        rm "${file}_sorted.tmp" "${file}_col4.tmp" "${file}_col6.tmp" "${file}_col7.tmp" "${file2}_sorted.tmp" "${file2}_col4.tmp" "${file2}_col6.tmp" "${file2}_col7.tmp"

    fi


  done < "$sample_list"


  # Add the header to the files
  echo -e $header > $numreads_file
  echo -e $header > $coverage_file
  echo -e $header > $meandepth_file
  echo -e $header > $numreads_file2
  echo -e $header > $coverage_file2
  echo -e $header > $meandepth_file2

  # Add the content of the temporary files to the final files
  cat "$temp_file_numreads" >> $numreads_file
  cat "$temp_file_coverage" >> $coverage_file
  cat "$temp_file_meandepth" >> $meandepth_file
  cat "$temp_file_numreads2" >> $numreads_file2
  cat "$temp_file_coverage2" >> $coverage_file2
  cat "$temp_file_meandepth2" >> $meandepth_file2


  # Remove the temporary files
  rm "${temp_file_numreads}" "${temp_file_coverage}" "${temp_file_meandepth}" "${temp_file_numreads2}" "${temp_file_coverage2}" "${temp_file_meandepth2}" 

fi


if [ "$step" == "4" ] || [ "$step" == "all" ]; then

  echo "Creating BED file with gene annotations for gene profiling"

  awk 'NR>1{print $2 "\t" $3 "\t" $4 "\t" $1 "\t0\t" ($5 == "f" ? "+" : ($5 == "r" ? "-" : $5))}' $final_annotation/nonredundant_plasmids-gene-calls.txt > "$final_annotation/nonredundant_plasmids-gene-calls.bed" 

fi 


if [ "$step" == "5" ] || [ "$step" == "all" ]; then

  # obtain files to run bedtools multicov
  filelist=$(ls $results_dir/*/*/*-profiling-sorted.bam | tr "\n" " ")

  echo "First indexing the bam files..."
  samtools index -M $filelist --threads $threads

fi 

if [ "$step" == "6" ] || [ "$step" == "all" ]; then

  # obtain files to run bedtools multicov
  filelist=$(ls $results_dir/*/*/*-profiling-sorted.bam | tr "\n" " ")

  tab_samplelist=$(echo $filelist | tr " " "\n" | sed 's/\.bam$//g' | sed 's/.*\///g' | sed 's/-profiling-sorted//g' | tr "\n" "\t")
  #create a header with the files that will be used (in the same order, since bedtools multicov does not handle it)
  echo -e "Plasmid_ID\tstart_pos\tend_pos\tgene_callers_id\tnothing\tdirection\t$tab_samplelist" > $results_dir/gene_profiling.tsv

  echo "Running bedtools multicov to extract the number of reads mapped to each gene region"
  bedtools multicov -bams $filelist -bed "$final_annotation/nonredundant_plasmids-gene-calls.bed" >> $results_dir/gene_profiling.tsv
  echo 'The output will be: plasmid	start_pos	end_pos	gene_callers_id	(empty)	direction	number_of_reads_mapped(one column per sample)'

fi 