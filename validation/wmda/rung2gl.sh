#!/bin/bash

export PERL_LWP_SSL_VERIFY_HOSTNAME=0
mkdir -p output
DONFILE_IN=../../graph_generator/data/wmda/don.txt
DONFILE_OUT=./output/don.gl.txt
./g2gl.pl -i $DONFILE_IN -o $DONFILE_OUT

PATFILE_IN=../../graph_generator/data/wmda/pat.txt
PATFILE_OUT=./output/pat.gl.txt
./g2gl.pl -i $PATFILE_IN -o $PATFILE_OUT

