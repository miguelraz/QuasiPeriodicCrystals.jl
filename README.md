# QuasiPeriodicCrystals

This repository contains a reimplementation of the algorithm published in [Efficient algorithm for simulating particles in true quasiperiodic environments by Alan Rodrigo Mendoza and Atahualpa S. Kraemer](https://arxiv.org/pdf/2111.08128.pdf) in collaboration with Miguel Raz Guzman Macedo.

To use this package, run 

```julia
using Pkg
Pkg.add("QuasiPeriodicCrystals")
```
after opening a [JuliaLang](https://julialang.org/downloads/) instance.

You should then be able to run
```julia
using QuasiPeriodicCrystals
sl = 1e6;
nsides = 5;
precision = BigFloat;
α = 0.0;
β = 15;
apoint = arb_Point(SL);
radiuscluster = 35.0;

MainClusterSites = main_Cluster(nsides, precision, α, β, apoint, radiuscluster)
```

This package is under heavy development and you should expect breaking changes unless you set version manually in your environment.
```
]add QuasiPeriodicCrystals#2.0.0
```