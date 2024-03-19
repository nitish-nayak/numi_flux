#!/bin/bash

# Make box around text @climagic
function box() { t="$1xxxx";c=${2:-=}; echo ${t//?/$c}; echo "$c $1 $c"; echo ${t//?/$c}; }

if [[ -z "$NUMIANA_DIR" ]]; then
    echo "numi_flux repo not setup. Run the setup script! Exiting!"
    exit 1
fi

# defaults
SEED=42
NJOBS=200
FILEPATH="/pnfs/uboone/persistent/users/bnayak/flux_files/uboone_geometrybugfix/me000z200i/run0/files"
OUTDIR="/pnfs/"${EXPERIMENT}"/scratch/users/"${USER}"/test_numiana"

usage() { echo "Usage: $0 [-p] [-s SEED] [-n NJOBS] [-i INPUT_FILELIST] [-o OUTPUT_DIR]" 1>&2; exit 1; }

while getopts 'ps:n:i:o:' args; do
    case "${args}" in
        p) echo "using local ppfx"
            local_ppfx=true ;;
        s)
            SEED="${OPTARG}" ;;
        n)
            NJOBS="${OPTARG}" ;;
        i)
            FILEPATH="${OPTARG}" ;;
        o)
            OUTDIR="${OPTARG}" ;;
        *)
            usage
            exit 1 ;;
    esac
done

box "Configuration : "
echo
echo "Seed : "$SEED
echo "Number of Jobs : "$NJOBS
echo "Input pnfs file path : "$FILEPATH
echo "Job Output directory : "$OUTDIR

TARDIR="/pnfs/"${EXPERIMENT}"/scratch/users/"${USER}"/grid_cache/numiana/"${RANDOM}
if [ ! -d "${TARDIR}" ]; then
    mkdir -p "${TARDIR}"
fi

cd "${NUMIANA_DIR}"
if [ ! -f "${TARDIR}"/local_numiflux.tar.gz ]; then
    echo "Tarring up numi_flux repo"
    tar -czf "${TARDIR}"/local_numiflux.tar.gz .
fi

if [ "$local_ppfx" = true ]; then
    if [[ -z "${PPFX_DIR}" ]]; then
        echo "PPFX_DIR is not set but requested local ppfx. Exiting!"
        exit 1
    fi
    if [ ! -f "${TARDIR}"/local_ppfx.tar.gz ]; then
        echo "Tarring up local ppfx"
        cd "${PPFX_DIR}"
        tar -czf "${TARDIR}"/local_ppfx.tar.gz .
        cd "${NUMIANA_DIR}"
    fi
fi

MACRO="jobsub/run_ppfxunivs.py"

if [ ! -d "${OUTDIR}" ]; then
    mkdir -p "${OUTDIR}"
fi

if [ ! -d "${PWD}"/jobsub/tmp ]; then
    mkdir -p "${PWD}"/jobsub/tmp
fi

echo "Making filelist from "${FILEPATH}" for "${NJOBS}
ls "${FILEPATH}"/*.root | head -n "${NJOBS}" | xargs pnfsToXRootD >> jobsub/tmp/input_files.txt

LOGFILE=${OUTDIR}"/numi_flux_log_\${PROCESS}_seed"${SEED}".log"

echo "Logfile is "${LOGFILE}
echo "Submitting jobs to "$OUTDIR
jobsub_submit --OS=SL7 \
              --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,ONSITE \
              --role=Analysis --memory=1900MB --expected-lifetime=48h --disk=10GB \
              -N "${NJOBS}" \
              -d NUMIANA "${OUTDIR}" \
              -G "${EXPERIMENT}" \
              -e MACRO="${MACRO}" \
              -e SEED="${SEED}" \
              -f "${TARDIR}"/local_numiflux.tar.gz \
              -f "${TARDIR}"/local_ppfx.tar.gz \
              -f dropbox://jobsub/tmp/input_files.txt \
              -f dropbox://"${MACRO}" \
              -L "${LOGFILE}" \
              file://jobsub/numiana_job.sh

echo "Cleaning up"
rm -r "${PWD}"/jobsub/tmp
