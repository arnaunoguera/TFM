#!/bin/bash


# Function to display help message
show_help() {
    echo "Usage: $0 -i <input_dir> -o <results_dir> -t <threads> -s <sample_list> [ -p <preset> -c ]"
    echo
    echo "Options:"
    echo "  -i  INPUT   Input directory (directory with the results of KNEADDATA, with one folder per sample inside)"
    echo "              Default: /mnt/synology/RAW_DATA/POP_STUDY/SYMLINKS_KNEADDATA"
    echo "  -o  OUTPUT  Output directory to store assemblies. Default: BASE_DIR/MEGAHIT/POP"
    echo "  -t  THREADS Number of threads (default: 1)"
    echo "  -s  SAMPLE LIST "
    echo "              File with a list of the samples to run this script for. "
    echo "              Default: BASE_DIR/samplelist.txt"
    echo "              (Unnecessary if -c is selected)"
    echo "  -p  PRESET  Use preset parameter for Megahit, i.e., meta-large (default: none)"
    echo "  -c          Perform coassembly with all samples in sample list, instead of one assembly per sample"
    echo "  -P  PSSWD   Password to run sudo commands"
    echo "  -h          Display this help message"
}

# Get options to run different parts of the strainphlan pipeline
while getopts i:o:t:s:p:c:P:h flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) results_dir=${OPTARG};;
        t) threads=${OPTARG};;
        s) sample_list=${OPTARG};;
        p) preset=${OPTARG};;
        c) coassembly=true;;
        P) passwd=${OPTARG};;
        h) show_help
           exit 0;;
        *) show_help
           exit 1;;
    esac
done

if [ -z "$threads" ]; then
  threads=1
fi

if [ -z "$input_dir" ]; then
  input_dir="/mnt/synology/RAW_DATA/POP_STUDY/SYMLINKS_KNEADDATA"
fi

if [ -z "$results_dir" ]; then
  results_dir="BASE_DIR/MEGAHIT/POP"
fi

if [ -z "$sample_list" ]; then
  sample_list="BASE_DIR/samplelist.txt"
fi


#####MAIN SCRIPT######
if [ -d $input_dir ]; then
  echo "Input directory is $input_dir"
else
  echo "ERROR: Couldn't find input directory. Exitting..."
  exit 1
fi

# Run MEGAHIT for each sample or as a coassembly depending on the flag
if [ "$coassembly" != true ]; then

  # Perform individual assemblies
  for sample in $(cat $sample_list);
  do

    echo "Starting process for $sample ..."
    if [ -d "$results_dir/${sample}/" ]; then
      echo "$results_dir/${sample}/ already exists, skipping to the next sample..."
      continue
    fi

    # Run MEGAHIT
    if [ -z "$preset" ]; then
      time megahit -1 "$input_dir/${sample}/${sample}_kneaddata_paired_1.fastq.gz" -2 "$input_dir/${sample}/${sample}_kneaddata_paired_2.fastq.gz" \
      -o "$results_dir/${sample}/" \
      --out-prefix "$sample" \
      --min-contig-len 1000 \
      --num-cpu-threads "$threads" 
    else 
      time megahit -1 "$input_dir/${sample}/${sample}_kneaddata_paired_1.fastq.gz" -2 "$input_dir/${sample}/${sample}_kneaddata_paired_2.fastq.gz" \
      -o "$results_dir/${sample}/" \
      --out-prefix "$sample" \
      --min-contig-len 1000 \
      --num-cpu-threads "$threads" \
      --presets "$preset"
    fi
  
  # Clear swap
  echo $passwd | sudo -S swapoff -a; sudo -S swapon -a

  done

else
  # Perform coassembly
  echo "Starting coassembly process..."
  if [ -d "$results_dir/COASSEMBLY/" ]; then
    echo "$results_dir/COASSEMBLY/ already exists, quitting..."
    exit 1
  fi

  #Get comma-separated lists of forwards and reverse reads
  R1s=`ls "$input_dir"/*/*_kneaddata_paired_1.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
  R2s=`ls "$input_dir"/*/*_kneaddata_paired_2.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`

  # Run MEGAHIT
  if [ -z "$preset" ]; then
    time megahit -1 $R1s -2 $R2s \
    -o "$results_dir/COASSEMBLY/" \
    --out-prefix "coassembly" \
    --min-contig-len 1000 \
    --num-cpu-threads "$threads" 
  else 
    time megahit -1 $R1s -2 $R2s \
    -o "$results_dir/COASSEMBLY/" \
    --out-prefix "coassembly" \
    --min-contig-len 1000 \
    --num-cpu-threads "$threads" \
    --presets "$preset"
  fi

fi
