#!/bin/sh
# run  comparison

OUTPUT_DIR=../../multi_race_impute/output
DATA_DIR=./data
date
python3 cmp_imp_results.py -o $DATA_DIR/pat.2015.txt.gz -n $OUTPUT_DIR/pat.gl.txt_out >$OUTPUT_DIR/pat.cmp_results.txt
python3 cmp_imp_results.py -o $DATA_DIR/don.2015.txt.gz -n $OUTPUT_DIR/don.gl.txt_out >$OUTPUT_DIR/don.cmp_results.txt

