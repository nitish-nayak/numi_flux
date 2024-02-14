This was originally a fork of Marco's original NuMIFlux [repo](https://github.com/marcodeltutto/NuMIFlux) for flugg files, combined with Krishan's [modifications](https://github.com/kvjmistry/NuMIFlux) for dk2nu files. However, it turns out that as I progressively keep making severe changes (including overhauling the code structure), this is worth detaching and maintaining separately. All credit to them and other predecessors for the main meat of the code.

- Various simplifications mainly
- Plan to improve some more areas to ensure easy comparisons
- Currently, `FluggFlux.cc` looks at FLUGG ntuples and `Dk2NuFlux.cc` looks at Dk2Nu ntuples
- Check `<flux_type>_example.py` in `python/` directory for example usage

## Instructions

After setting up `uboonecode` or optionally only the dependencies (see below)
```
source setup_numiana.sh
make all
```
If you just want to work with individual `<flux_type>`s where `<flux_type>` is either `dk2nu` or `flugg`, one can do
```
make base
make <flux_type>
```
To start afresh, as usual one can do
```
make clean
```

### Dependencies
- `flugg` just needs `root` basically
- `dk2nu` pulls in additionally `ppfx`, `boost`, `dk2nu`

### File Locations
- `flugg` : originally available on bluearc
    - FHC : `/nusoft/data/flux/blackbird-numix/flugg_mn000z200i_rp11_lowth_pnut_f11f093bbird_target/root`
    - RHC : `/nusoft/data/flux/blackbird-numix/flugg_mn000z-200i_rp11_lowth_pnut_f11f093bbird_target/root`
    - Also copied to `/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files` and `/pnfs/uboone/persistent/users/bnayak/flux_files/flugg_files` (dcache)
    - See [wiki](https://cdcvs.fnal.gov/redmine/projects/numi-beam-sim/wiki/Locations_of_shared_files) for more information

- `dk2nu` : available on stashcache
    - FHC : `/cvmfs/uboone.osgstorage.org/stash/uboonebeam/numi_dk2nu_zero_threshold/FHC`
    - RHC : `/cvmfs/uboone.osgstorage.org/stash/uboonebeam/numi_dk2nu_zero_threshold/RHC`

## Planned Changes

- ~Common output interface, currently its replicated across both `<flux_type>` classes~
- Parallel Processing
    - Atleast for FLUGG, this should be easy
    - Can do for Dk2Nu as well but PPFX has a bunch of static instances that is likely not thread safe
