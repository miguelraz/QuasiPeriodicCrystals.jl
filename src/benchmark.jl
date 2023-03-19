using LinearAlgebra
using MutableArithmetics
using Random
using StaticArrays
using Test
#include("voronoi.jl")
include("q4.jl")
#include("quasi2.jl")
#include("Quasicrystals.jl")
#include("qset.jl")
#using QuasiPeriodicCrystals
Random.seed!(42)
SL = 1e6
NSides = 7
Precision = BigFloat;
α = 0.0
# \beta = 34 # also works!
β = 40;
APoint = arb_Point(SL)
RadiusCluster = 120.0;

# 3.66s, 52M allocs, 2.1Gb, 33% GC time
# length(MainClusterSites) == 56625
@show Precision
@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster);
@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster);
@show length(MainClusterSites)
@assert length(MainClusterSites) == 56625 "BigFloat FAIL"
Precision = Float64
@show Precision
@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster);
@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster);
@show length(MainClusterSites)
@assert length(MainClusterSites) == 152003 "Float64 FAIL"
#Precision = Float32
#@show Precision
#@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster);
#@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster);
#@show length(MainClusterSites)
#@assert length(MainClusterSites) == 152003 "Float64 FAIL"
