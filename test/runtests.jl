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



inputs = [[parse(Float64, x) for x in eachline(s)] for s in APoint]
@test length(inputs) == 10

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

outputs = [[parse.(BigFloat, split(s, ',')) for s in eachline(x)] for x in Sites]

@test length(outputs) == 10

@testset "QuasiPeriodicCrystals.jl" begin
    # Write your tests here.
    @test dot([1, 2], [3, 4]) == 11
    @test dot_Product([1, 2], [3, 4]) == 11
    @test orthogonal_Vec([1, 2]) == (2, -1)
    for _ in 1:100
        @test hypot.(arb_Point(3.0)...) <= 3√2
    end
end


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
RTests_inputs = [
    "Sites_Test_N5.csv"
    "Sites_Test_N8.csv"
    "Sites_Test_N23.csv"
]

RTests_outputs = [[parse.(BigFloat, split(s, ',')) for s in eachline(x)] for x in RTests_inputs]
@test length(RTests_outputs) == 3
@test length(RTests_outputs[1]) == 4739
@test length(RTests_outputs[2]) == 4651
@test length(RTests_outputs[3]) == 4888

