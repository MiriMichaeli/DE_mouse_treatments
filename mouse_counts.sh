#!/bin/bash

# https://www.gencodegenes.org/mouse/release_M10.html download the GFF3 from "Comprehensive gene annotation"
# save both the gff and the input.sam on your computer (you will need to insert the full paths as input to this script)
# you need the HTSeq tool installed on your computer and added to PATH
# how to run: sh mouse_counts.sh <input_sam_path> <gff3_path>

# run HTSeq to count reads per genes:
htseq-count --stranded=no --additional-attr=gene_name $1 $2 > counts_out.txt
# remove the last lines of counters:
grep -v "__" counts_out.txt > tmpfile && mv tmpfile counts_out.txt 
# prepare the output fille with gene name and counts columns and filter to genes with > 10 reads:
cat counts_out.txt | awk '{print $2"\t"$3}' | awk '{ if ($2 > 10) { print } }' > counts_filtered.txt
