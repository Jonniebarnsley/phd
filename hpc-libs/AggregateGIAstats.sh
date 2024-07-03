#!/bin/bash

# script to concatenate all individual statsfiles into a single summary file

usage() { echo "Usage: $0 <GIAstats_path> <destination_filepath>" 1>&2; exit 1; }

# usage clause
if [ "$#" -ne 2 ]; then
    usage
fi

GIAstats="$1"
outfile="$2"

# if the summary already contains data up to 9990 years, stop
if [ -f $outfile ] && tail $outfile | grep -q 'time = 9.990'; then
    exit
fi

# otherwise, wipe the summary file and copy all statsfiles into summary
> "$outfile"
for statsfile in $GIAstats/*.GIAstats; do
    grep 'time = ' "$statsfile" >> "$outfile"
done
 
