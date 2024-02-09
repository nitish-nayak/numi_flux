This was originally a fork of Marco's original NuMIFlux [repo](https://github.com/marcodeltutto/NuMIFlux) for flugg files, combined with Krishan's [modifications](https://github.com/kvjmistry/NuMIFlux) for dk2nu files. However, it turns out that as I progressively keep making severe changes (including overhauling the code structure), this is worth detaching. All credit to them and other predecessors for the main meat of the code.

- Various simplifications mainly
- Plan to improve some more areas to ensure easy comparisons
- Currently, `FluggFlux.cc` looks at FLUGG ntuples and `Dk2NuFlux.cc` looks at Dk2Nu ntuples
- Check `<flux_type>_example.py` in `python/` directory for example usage

## Instructions

```
source setup_uboone.sh
make all
```
If you just want to work with individual `<flux_type>`s where `<flux_type>` is either `dk2nu` or `flugg`, one can do
```
make <flux_type>
```
To start afresh, as usual one can do
```
make clean
```

### Dependencies
- `flugg` just needs `root` basically
- `dk2nu` pulls in additionally `ppfx`, `boost`, `dk2nu`

## Planned Changes

- ~Common output interface, currently its replicated across both `<flux_type>` classes~
- Parallel Processing
    - Atleast for FLUGG, this should be easy
    - Can do for Dk2Nu as well but PPFX has a bunch of static instances that is likely not thread safe
