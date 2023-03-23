using QuasiPeriodicCrystals
using LinearAlgebra
using StaticArrays
using Test

@testset "QuasiPeriodicCrystals.jl" begin
    # Write your tests here.
    @test dot([1, 2], [3, 4]) == 11
    @test dot_Product([1, 2], [3, 4]) == 11
    @test orthogonal_Vec([1, 2]) == (2, -1)
    for _ in 1:100
        @test hypot.(arb_Point(3.0)...) <= 3√2
    end
end
