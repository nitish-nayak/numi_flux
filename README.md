# NuMIFlux
Code to calculate the NuMI Flux at MicroBooNE


Run with:

```
source SetupNuMIFlux.sh

root -l
> gSystem->Load("FluggNtuple/FluxNtuple_C.so");
> gSystem->Load("NuMIFlux_cc.so");
> NuMIFlux f;
> f.CalculateFlux();
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
