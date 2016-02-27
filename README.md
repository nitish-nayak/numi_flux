# NuMIFlux
Code to calculate the NuMI Flux at MicroBooNE

Setup ROOT:

```
source /nusoft/app/externals/setup
setup root   v5_34_18a   -q e5:debug:nu
```

Run with:

```
root -l
> gSystem->Load("NuMIFlux_cc.so");
> NuMIFlux f;
> f.CalculateFlux();
```
