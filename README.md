# NuMIFlux
Code to calculate the NuMI Flux at MicroBooNE.

Wiki page: https://cdcvs.fnal.gov/redmine/projects/ubooneoffline/wiki/NuMI_Flux_Histograms.


## First Usage

From the MicroBooNE gpvm machines:

```
git clone git@github.com:marcodeltutto/NuMIFlux.git
cd NuMIFlux/
source SetupNuMIFlux.sh

cd FluggNtuple/
root -l 
> .L FluxNtuple.C+
> .q

cd ../
root -l
> gSystem->Load("FluggNtuple/FluxNtuple_C.so");
> .L NuMIFlux.cc+
> .q

```

## General Instructions

Run with:

```
source SetupNuMIFlux.sh

root -l
> gSystem->Load("FluggNtuple/FluxNtuple_C.so");
> gSystem->Load("NuMIFlux_cc.so");
> NuMIFlux f("generic_*_to_flugg.root");
> f.CalculateFlux();
```

or use the python script:

```
source SetupNuMIFlux.sh
python RunNuMIFlux.py
```

Compile with:

```
root -l
> gSystem->Load("FluggNtuple/FluxNtuple_C.so");
> .L NuMIFlux.cc+
```


Generate new FluxNtuple_C.so:

```
cd FluggNtuple
root -l
.L FluxNtuple.C+
```
