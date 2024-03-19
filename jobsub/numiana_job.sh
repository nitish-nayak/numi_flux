#!/bin/bash

# Make box around text @climagic
function box() { t="$1xxxx";c=${2:-=}; echo ${t//?/$c}; echo "$c $1 $c"; echo ${t//?/$c}; }

function setup_localppfx {
    cd $CONDOR_DIR_INPUT
    mkdir -p "$CONDOR_DIR_INPUT"/local_ppfx
    mv local_ppfx.tar.gz local_ppfx/.
    cd "$CONDOR_DIR_INPUT"/local_ppfx

    echo "Untarring local ppfx"
    tar xvfz local_ppfx.tar.gz -C ./ > /dev/null

    echo "Setting up local PPFX"
    export PPFX_DIR=${PWD}
    source setup_ppfx_grid.sh

    cd $CONDOR_DIR_INPUT
}

box "At CONDOR_DIR_INPUT"
cd $CONDOR_DIR_INPUT
pwd
ls

echo
box "Untarring numi_flux"
tar xvfz local_numiflux.tar.gz -C ./ > /dev/null

echo
box "Setting up various dependencies"
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup uboonecode v08_00_00_76 -q e17:prof

if [ -f local_ppfx.tar.gz ]; then
    echo
    box "Setting up local PPFX"
    setup_localppfx
    echo "PPFX_DIR : "$PPFX_DIR
    if [[ -z "${PPFX_DIR}" ]]; then
        echo "Error in setting up PPFX_DIR"
        exit
    fi
fi

echo
box "Setting up local numi_flux repo and building dictionaries"
source setup_numiana_grid.sh
make all

if [[ ! -f "${MACRO}" ]]; then
    echo "No python macro found! Exiting.."
    exit
fi

echo
FILELIST="input_files.txt"
box "The input file list is"
cat "$FILELIST"
echo "Splitting filelist into chunks of 1 file"
cat "${FILELIST}" | split -l 1 -a 3 -d
job_filelist="x"`printf "%03d" "${PROCESS}"`

box "Filelist for job is : "
cat "$job_filelist"

echo
box "Running python macro for process : "${PROCESS}
python "${MACRO}" "${job_filelist}" "${SEED}"

echo
box "Copying output file"
mv numi_flux_output.root "${CONDOR_DIR_NUMIANA}"/numi_flux_output_"${PROCESS}"_seed"${SEED}".root
