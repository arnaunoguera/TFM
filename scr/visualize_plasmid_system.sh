#!/bin/bash


# Script to: geta  pdf visualization of the gene content and identity of a group of contigs / plasmid system 
#######WARNING######
# This sctipt neds the conda environment mobmess activated

# Function to display help message
show_help() {
    echo "Usage: $0 -p <PLASMIDS_DOIR> -m <MOBMESS_DIR> -t <threads> [ -s <system> or -c <CONTIGS> ] [ -a -b <BLOCK_HEIGHT> -h ]"
    echo
    echo "Options:"
    echo "  -p  PLASMIDS_DOIR   Input dir (with redundant plasmids fasta file)"
    echo "                      Default: BASE_DIR/PLASMIDS"
    echo "  -g  GENE_ANNOTATIONS DIR with gene annotations (COGs and Pfams)"
    echo "                      Default: BASE_DIR/ANVIO/FINAL_PLASMIDS"
    echo "  -m  MOBMESS DIR     Directory with results of MobMess. Default: BASE_DIR/MOBMESS"
    echo "  -s  SYSTEM          PLASMID SYSTEM to visualize (ex: PS4599) (incompatible with -c)"
    echo "  -c  CONTIGS         Contigs to visualize (ex: HCPOP_420_000000011069,HCPOP_420_0000000110790) (incompatible with -s)"
    echo "  -a  ALL             Visualize all plasmids in a system (instead of just the represenantives of each cluster) (default: false)"
    echo "  -b  BLOCK_HEIGHT    Height of the blocks in the visualization (default: 0.5)"
    echo "  -t  THREADS         Number of threads to use (default: 1)"
    echo "  -h          Display this help message"
}

# Initialize variables
all=0

# Get options to run different parts of the strainphlan pipeline
while getopts p:m:g:s:b:c:t:ah flag
do
    case "${flag}" in
        p) plasmid_dir=${OPTARG};;
        m) mobmess=${OPTARG};;
        s) system=${OPTARG};;
        c) contigs=${OPTARG};;
        t) threads=${OPTARG};;
        a) all=1;;
        b) block_height=${OPTARG};;
        g) gene_annotations=${OPTARG};;
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

if [ -z "$mobmess" ]; then
  mobmess="BASE_DIR/MOBMESS"
fi

if [ -z "$system" ] && [ -z "$contigs" ]; then
    echo "Error: Either -s or -c must be provided."
    show_help
    exit 1
elif [ -n "$system" ] && [ -n "$contigs" ]; then
    echo "Error: Both -s and -c cannot be provided at the same time."
    show_help
    exit 1
fi

if [ -z "$threads" ]; then
  threads="1"
fi

if [ -z "$block_height" ]; then
  block_height="0.5"
fi

if [ -z "$gene_annotations" ]; then
  gene_annotations="BASE_DIR/ANVIO/FINAL_PLASMIDS"
fi


#### MAIN SCRIPT 

# If system is provided, extract list of contigs in the system
if [ -n "$system" ]; then

    if [ $all -eq 1 ]; then
        # If all contigs should be visualized, just get all backbone ($4) and compound ($8) plasmids
        contigs=$(awk -v sys="$system" '$1 == sys {print $4 "|" $8 }'  "$mobmess/plasmids-mobmess_systems.txt" | tr "|" ",")
    else 
        # create a tmp folder if it doesn't already exist in MOBMESS
        if [ ! -d "$mobmess/tmp" ]; then
            mkdir -p "$mobmess/VISUALIZATION/tmp"
        fi
        
        # get the clusters first: backbone cluster ($2) and compound clusters ($5), to then obtain the representative of each cluster
        awk -v sys="$system" '$1 == sys {print $2 "|" $5 }'  "$mobmess/plasmids-mobmess_systems.txt" | tr "|" "\n" > "$mobmess/VISUALIZATION/tmp/clusters.txt"
        # get the representative of each cluster
        contigs=$(awk 'FNR==NR { ids[$1]; next } $2 in ids {print $5}' "$mobmess/VISUALIZATION/tmp/clusters.txt" "$mobmess/plasmids-mobmess_clusters.txt" | sort | uniq | paste -d "," -s)
    fi

else
    # If contigs are provided just use them
    # But save the system variable for the file name of the pdf as: contigs_currentdate
    system=$(date +"contigs-%Y-%m-%d")

fi

echo "Visualizing system $system with contigs: $contigs"

# Extract only the relevant sequences from the fasta file 
# First get the list of contigs and then use seqtk to extract them
echo -e $contigs | tr "," "\n" > "$mobmess/VISUALIZATION/tmp/contigs.txt"
# Then extract the fasta 
grep --no-group-separator -f $mobmess/VISUALIZATION/tmp/contigs.txt $plasmid_dir/predicted_plasmids_redundant.fa -A 1 > "$mobmess/VISUALIZATION/tmp/contigs.fa"



# Run mobmess itself
mobmess visualize --sequences "$mobmess/VISUALIZATION/tmp/contigs.fa" --annotations $gene_annotations/nonredundant_plasmids-cogs-and-pfams.txt \
 --contigs $contigs --threads $threads --output "$mobmess/VISUALIZATION/${system}.pdf" --gene-calls $gene_annotations/nonredundant_plasmids-gene-calls.txt \
 --align-blocks-height $block_height

echo "Results in $mobmess/VISUALIZATION/${system}.pdf"