#!/bin/bash

# files
infile=$1
outfile=$2

# input values
flowcell="2436449055"
lane="1"

awk -F'\t' -v fc_val="$flowcell" -v lane_val="$lane" '
  NR==1 {
    for(i=1;i<=NF;i++) {
      if($i=="Flowcell") fs_col=i;
      if($i=="Lane") lane_col=i;
    }
    # Verify we found both columns
    if(!fs_col || !lane_col) {
      print "Error: Could not find required column headers" > "/dev/stderr";
      exit 1;
    }
    print;
    next
  }
  $fs_col == fc_val && $lane_col == lane_val
' "$infile" > "$outfile"
