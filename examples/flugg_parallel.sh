#!/bin/bash

OUTDIR="."

usage() { echo "Usage: $0 [-n NTHREADS] [-o OUTDIR] [INPUT_FILELIST]" 1>&2; exit 1; }

# Make box around text @climagic
function box() { t="$1xxxx";c=${2:-=}; echo ${t//?/$c}; echo "$c $1 $c"; echo ${t//?/$c}; }

while getopts ':n:o:' args; do
    case "${args}" in
        n)
            NTHREADS="${OPTARG}" ;;
        o)
            OUTDIR="${OPTARG}" ;;
        *|h)
            usage
            exit 1 ;;
    esac
done

shift $((OPTIND - 1))

if [ "${NTHREADS}" -eq "${NTHREADS}" ] 2>/dev/null; then
    INPUT_FILELIST="${@}"
    box "Configuration : "
    echo "NTHREADS : "${NTHREADS}
    echo "OUTDIR : "${OUTDIR}
else
    echo "Didn't get numeric input for NTHREADS. Exiting!"
    usage
    exit
fi

tr ' ' '\n' <<< "${INPUT_FILELIST}" >> tmp_filelist.txt
NFILES=`cat tmp_filelist.txt | wc -l`
if [ "$NFILES" -lt 1 ]; then
    box "Zero input files, malformed input probably. Exiting!"
    rm tmp_filelist.txt
    exit
fi

cat tmp_filelist.txt | split -n r/"${NTHREADS}" -a 2 -d

cat << EOF > tmp_flugg_parallel.py
import os
import sys
import re
import ROOT

ROOT.gSystem.Load("\$NUMIANA_DIR/lib/FluggFlux_cc.so")
ROOT.gROOT.SetBatch(True)

filelist = sys.argv[1]
outfile = sys.argv[2]

print("Reading from : %s" % filelist)
print("Output file is : %s" % outfile)

f = ROOT.FluggFlux(True, filelist, outfile)
f.CalculateFlux()
EOF

for thread in `seq 0 $((${NTHREADS}-1))`; do
    box "Starting Thread : "$thread
    python tmp_flugg_parallel.py x0"${thread}" "${OUTDIR}"/flugg_output_numiana_part"${thread}".root &
    sleep 5
done

wait

box "Cleaning up"
rm tmp_filelist.txt
rm tmp_flugg_parallel.py
rm x0*
