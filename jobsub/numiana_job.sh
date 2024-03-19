#!/bin/bash

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


echo "At CONDOR_DIR_INPUT"
cd $CONDOR_DIR_INPUT
pwd
ls

echo "Untarring numi_flux"
tar xvfz local_numiflux.tar.gz -C ./ > /dev/null

echo "Setting up various dependencies"
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup uboonecode v08_00_00_76 -q e17:prof

if [ -f local_ppfx.tar.gz ]; then
    echo "Setting up local PPFX"
    setup_localppfx
    echo "PPFX_DIR : "$PPFX_DIR
    if [[ -z "${PPFX_DIR}" ]]; then
        echo "Error in setting up PPFX_DIR"
        exit
    fi
fi

echo "Setting up local numi_flux repo"
source setup_numiana_grid.sh

if [[ ! -f "${MACRO}" ]]; then
    echo "No python macro found! Exiting.."
    exit
fi

echo "Running python macro"
python "${MACRO}" "${FILELIST}"

echo "Copying output file"
mv numi_flux_output.root "${CONDOR_DIR_NUMIANA}"/numi_flux_output_"${PROCESS}".root
