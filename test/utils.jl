using ExactOptimalTransport

using LinearAlgebra
using Random
using Test
using Distributions

Random.seed!(100)

@testset "utils.jl" begin
    @testset "FiniteDiscreteMeasure" begin
        @testset "Univariate Finite Discrete Measure" begin
            n = 100
            m = 80
            μsupp = rand(n)
            νsupp = rand(m)
            μprobs = normalize!(rand(n), 1)

            μ = ExactOptimalTransport.discretemeasure(μsupp, μprobs)
            ν = ExactOptimalTransport.discretemeasure(νsupp)
            # check if it vectors are indeed probabilities
            @test isprobvec(μ.p)
            @test isprobvec(probs(μ))
            @test ν.p == ones(m) ./ m
            @test probs(ν) == ones(m) ./ m

            # check if it assigns to DiscreteNonParametric when Vector/Matrix is 1D
            @test μ isa DiscreteNonParametric
            @test ν isa DiscreteNonParametric

            # check if support is correctly assinged
            @test sort(μsupp) == μ.support
            @test sort(μsupp) == support(μ)
            @test sort(vec(νsupp)) == ν.support
            @test sort(vec(νsupp)) == support(ν)
        end
        @testset "Multivariate Finite Discrete Measure" begin
            n = 10
            m = 3
            μsupp = [rand(m) for i in 1:n]
            νsupp = [rand(m) for i in 1:n]
            μprobs = normalize!(rand(n), 1)
            μ = ExactOptimalTransport.discretemeasure(μsupp, μprobs)
            ν = ExactOptimalTransport.discretemeasure(νsupp)
            # check if it vectors are indeed probabilities
            @test isprobvec(μ.p)
            @test isprobvec(probs(μ))
            @test ν.p == ones(n) ./ n
            @test probs(ν) == ones(n) ./ n

            # check if support is correctly assinged
            @test μsupp == μ.support
            @test μsupp == support(μ)
            @test νsupp == ν.support
            @test νsupp == support(ν)
        end
    end
end
