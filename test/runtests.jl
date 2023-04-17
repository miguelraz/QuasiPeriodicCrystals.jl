using QuasiPeriodicCrystals
using LinearAlgebra
using StaticArrays
using Test

APoint = [
    "APoint_Test_1.csv",
    "APoint_Test_2.csv",
    "APoint_Test_3.csv",
    "APoint_Test_4.csv",
    "APoint_Test_5.csv",
    "APoint_Test_6.csv",
    "APoint_Test_7.csv",
    "APoint_Test_8.csv",
    "APoint_Test_9.csv",
    "APoint_Test_10.csv"
]

APoint_inputs = [[parse(Float64, x) for x in eachline(s)] for s in APoint]
@test length(APoint_inputs) == 10
@test all(==(2), length.(APoint_inputs))

Sites = [
    "Sites_Test_1.csv"
    "Sites_Test_2.csv"
    "Sites_Test_3.csv"
    "Sites_Test_4.csv"
    "Sites_Test_5.csv"
    "Sites_Test_6.csv"
    "Sites_Test_7.csv"
    "Sites_Test_8.csv"
    "Sites_Test_9.csv"
    "Sites_Test_10.csv"
]

Sites_outputs = [[parse.(BigFloat, split(s, ',')) for s in eachline(x)] for x in Sites]

@test length(Sites_outputs) == 10
@test length(Sites_outputs[1]) == 4897
@test length(Sites_outputs[2]) == 4880
@test length(Sites_outputs[3]) == 4917
@test length(Sites_outputs[4]) == 4897
@test length(Sites_outputs[5]) == 4880
@test length(Sites_outputs[6]) == 4902
@test length(Sites_outputs[7]) == 4894
@test length(Sites_outputs[8]) == 4913
@test length(Sites_outputs[9]) == 4899
@test length(Sites_outputs[10]) == 4889

#@testset "QuasiPeriodicCrystals.jl" begin
# Write your tests here.
@test dot([1, 2], [3, 4]) == 11
#@test dot_Product([1, 2], [3, 4]) == 11
#@test orthogonal_Vec([1, 2]) == (2, -1)
#end

@test 1 == 1

#=
# Tests de Alan

# Test 1
import Random
Random.seed!(1234)

SL = 1e6;
NSides = 5;
Precision = BigFloat;
α = 0.0;
β = 15;
APoint = arb_Point(SL);
RadiusCluster = 35.0;

@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster);

writedlm("Sites_Test_N$(NSides).csv", MainClusterSites, ',');

## Test 2
import Random
Random.seed!(1234)

SL = 1e6;
NSides = 8;
Precision = BigFloat;
α = 0.0;
β = 10;
APoint = arb_Point(SL);
RadiusCluster = 35.0;

@time MainClusterSites = main_Cluster(NSides, Precision, α, β, APoint, RadiusCluster);

writedlm("Sites_Test_N$(NSides).csv", MainClusterSites, ',');

# Test 3

=#
RTests_outputs = [
    "Sites_Test_N5.csv"
    "Sites_Test_N8.csv"
    "Sites_Test_N23.csv"
]

RTests_outputs = [[parse.(BigFloat, split(s, ',')) for s in eachline(x)] for x in RTests_outputs]
@test length(RTests_outputs) == 3
@test length(RTests_outputs[1]) == 4739
@test length(RTests_outputs[2]) == 4651
@test length(RTests_outputs[3]) == 4888

