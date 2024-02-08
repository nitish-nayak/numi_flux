Fork of original NuMIFlux [repo](https://github.com/marcodeltutto/NuMIFlux), combined with Krishan's [modifications](https://github.com/kvjmistry/NuMIFlux) for dk2nu files

- Various simplifications as well
- Plan to improve some areas to ensure easy comparisons
- Currently, FluggFlux.cc looks at FLUGG ntuples and Dk2NuFlux.cc looks at Dk2Nu ntuples
- Check RunFluggFlux.py and RunDk2NuFlux.py in the `<flux_type>` (see below) directory for usage

## Instructions

```
source setup_uboone.sh
cd <flux_type>
make <flux_type>
```
where `<flux_type>` is either `dk2nu` or `flugg`.

## Planned Changes

- Common output interface, currently its replicated across both `<flux_type>` classes
- Parallel Processing
    - Atleast for FLUGG, this should be easy
    - Can do for Dk2Nu as well but PPFX has a bunch of static instances that is likely not thread safe
