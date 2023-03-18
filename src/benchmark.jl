include("voronoi.jl")
include("Quasicrystals.jl")
#using QuasiPeriodicCrystals
using LinearAlgebra
using Random
Random.seed!(42)
SL = 1e6
NSides = 7
Precision = BigFloat;
α = 0.0
β = 40;
APoint = arb_Point(SL)
RadiusCluster = 120.0;
# 4.25s, 52M allocs, 2.1Gb, 33% GC time
@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster)