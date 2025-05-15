#!/bin/bash

#######WARNING######
# This sctipt neds the conda environment anvio8 activated

###### DEFINE VARIABLES TO RUN ANVI'O ######

data_dir="BASE_DIR/ANVIO"
threads="30"

#####MAIN SCRIPT######

databases=$(find $data_dir -name "*db" -type f)


anvi-display-contigs-stats $databases \
    --report-as-text \
    -o "BASE_DIR/ANVIO/contigs_stats.txt" 