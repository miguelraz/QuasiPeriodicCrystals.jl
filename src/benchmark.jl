using LinearAlgebra
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
β = 40;
APoint = arb_Point(SL)
RadiusCluster = 120.0;

# 3.66s, 52M allocs, 2.1Gb, 33% GC time
# length(MainClusterSites) == 56625
@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster);
@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster);
@show length(MainClusterSites)
@assert length(MainClusterSites) == 56625 "BigFloat does not pass"
@show "BigFloat"
Precision = Float64
@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster);
@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster);
@assert length(MainClusterSites) == 152003 "BigFloat passes"
