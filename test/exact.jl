using ExactOptimalTransport

using Distances
using FillArrays
using GLPK
using PythonOT: PythonOT
using Tulip
using MathOptInterface
using Distributions
using HCubature

using LinearAlgebra
using Random
using SparseArrays

const MOI = MathOptInterface
const POT = PythonOT

Random.seed!(100)

@testset "exact.jl" begin
    @testset "Earth-Movers Distance" begin
        M = 200
        N = 250
        μ = normalize!(rand(M), 1)
        ν = normalize!(rand(N), 1)

        @testset "example" begin
            # create random cost matrix
            C = pairwise(SqEuclidean(), rand(1, M), rand(1, N); dims=2)

            # compute optimal transport map and cost with POT
            pot_P = POT.emd(μ, ν, C)
            pot_cost = POT.emd2(μ, ν, C)

            # compute optimal transport map and cost with Tulip and GLPK
            for T in (Tulip.Optimizer, GLPK.Optimizer)
                lp = T()
                P = emd(μ, ν, C, lp)
                @test size(C) == size(P)
                @test MOI.get(lp, MOI.TerminationStatus()) == MOI.OPTIMAL
                @test maximum(abs, P .- pot_P) < 1e-2

                lp = T()
                cost = emd2(μ, ν, C, lp)
                @test dot(C, P) ≈ cost atol = 1e-5
                @test MOI.get(lp, MOI.TerminationStatus()) == MOI.OPTIMAL
                @test cost ≈ pot_cost atol = 1e-5
            end
        end

        @testset "pre-computed plan" begin
            # create random cost matrix
            C = pairwise(SqEuclidean(), rand(1, M), rand(1, N); dims=2)

            # compute optimal transport map with Tulip and GLPK
            for T in (Tulip.Optimizer, GLPK.Optimizer)
                P = emd(μ, ν, C, T())

                # do not use μ and ν to ensure that provided map is used
                cost = emd2(similar(μ), similar(ν), C, T(); plan=P)
                @test cost ≈ emd2(μ, ν, C, T())
            end
        end

        # https://github.com/JuliaOptimalTransport/OptimalTransport.jl/issues/71
        @testset "cost matrix with integers" begin
            C = pairwise(SqEuclidean(), rand(1:10, 1, M), rand(1:10, 1, N); dims=2)
            emd2(μ, ν, C, Tulip.Optimizer())
        end
    end

    @testset "1D Optimal Transport for Convex Cost" begin
        @testset "continuous distributions" begin
            # two normal distributions (has analytical solution)
            μ = Normal(randn(), 1 + rand())
            ν = Normal(randn(), 1 + rand())

            # compute OT plan
            γ = ot_plan(sqeuclidean, μ, ν)
            for x in randn(10)
                @test γ(x) ≈ invlogcdf(ν, logcdf(μ, x))
            end

            # compute OT cost
            c = ot_cost(sqeuclidean, μ, ν)
            @test c ≈ (mean(μ) - mean(ν))^2 + (std(μ) - std(ν))^2

            # do not use ν to ensure that the provided plan is used
            @test ot_cost(sqeuclidean, μ, Normal(randn(), rand()); plan=γ) ≈ c
        end

        @testset "semidiscrete case" begin
            μ = Normal(randn(), rand())
            νprobs = normalize!(rand(30), 1)
            ν = Categorical(νprobs)

            # compute OT plan
            γ = ot_plan(euclidean, μ, ν)
            for x in randn(10)
                @test γ(x) ≈ invlogcdf(ν, logcdf(μ, x))
            end

            # compute OT cost, without and with provided plan
            # do not use ν in the second case to ensure that the provided plan is used
            c = ot_cost(euclidean, μ, ν)
            @test ot_cost(euclidean, μ, Categorical(reverse(νprobs)); plan=γ) ≈ c

            # check that OT cost is consistent with OT cost of a discretization
            m = 500
            xs = rand(μ, m)
            μdiscrete = fill(1 / m, m)
            C = pairwise(Euclidean(), xs', (1:length(νprobs))'; dims=2)
            for optimizer in (Tulip.Optimizer(), GLPK.Optimizer())
                c2 = emd2(μdiscrete, νprobs, C, optimizer)
                @test c2 ≈ c rtol = 1e-1
            end
        end

        @testset "discrete case" begin
            # different random sources and target marginals:
            # non-uniform + different size, uniform + different size, uniform + equal size
            for (μ, ν) in (
                (
                    DiscreteNonParametric(randn(30), normalize!(rand(30), 1)),
                    DiscreteNonParametric(randn(50), normalize!(rand(50), 1)),
                ),
                (
                    DiscreteNonParametric(randn(30), Fill(1 / 30, 30)),
                    DiscreteNonParametric(randn(50), Fill(1 / 50, 50)),
                ),
                (
                    DiscreteNonParametric(randn(30), Fill(1 / 30, 30)),
                    DiscreteNonParametric(randn(30), Fill(1 / 30, 30)),
                ),
            )
                # extract support, probabilities, and "size"
                μsupport = support(μ)
                μprobs = probs(μ)
                m = length(μprobs)

                νsupport = support(ν)
                νprobs = probs(ν)
                n = length(νprobs)

                # compute OT plan
                γ = @inferred(ot_plan(euclidean, μ, ν))
                @test γ isa SparseMatrixCSC
                @test size(γ) == (m, n)
                @test vec(sum(γ; dims=2)) ≈ μ.p
                @test vec(sum(γ; dims=1)) ≈ ν.p

                # consistency checks
                I, J, W = findnz(γ)
                @test all(w > zero(w) for w in W)
                @test sum(W) ≈ 1
                @test sort(unique(I)) == 1:m
                @test sort(unique(J)) == 1:n
                @test sort(I .+ J) == if μprobs isa Fill && νprobs isa Fill && m == n
                    # Optimized version for special case (discrete uniform + equal size)
                    2:2:(m + n)
                else
                    # Generic case (not optimized)
                    2:(m + n)
                end

                # compute OT cost
                c = @inferred(ot_cost(euclidean, μ, ν))

                # compare with computation with explicit cost matrix
                # DiscreteNonParametric sorts the support automatically, here we have to sort
                # manually
                C = pairwise(Euclidean(), μsupport', νsupport'; dims=2)
                for optimizer in (Tulip.Optimizer(), GLPK.Optimizer())
                    c2 = emd2(μprobs, νprobs, C, optimizer)
                    @test c2 ≈ c rtol = 1e-5
                end

                # compare with POT
                # disabled currently since https://github.com/PythonOT/POT/issues/169 causes bounds
                # error
                # @test γ ≈ POT.emd_1d(μ.support, ν.support; a=μ.p, b=μ.p, metric="euclidean")
                # @test c ≈ POT.emd2_1d(μ.support, ν.support; a=μ.p, b=μ.p, metric="euclidean")

                # do not use the probabilities of μ and ν to ensure that the provided plan is
                # used
                μ2 = DiscreteNonParametric(μsupport, reverse(μprobs))
                ν2 = DiscreteNonParametric(νsupport, reverse(νprobs))
                c2 = @inferred(ot_cost(euclidean, μ2, ν2; plan=γ))
                @test c2 ≈ c
                c2 = @inferred(ot_cost(euclidean, μ2, ν2; plan=Matrix(γ)))
                @test c2 ≈ c
            end
        end
    end

    @testset "Multivariate Gaussians" begin
        @testset "translation with constant covariance" begin
            m = randn(100)
            τ = rand(100)
            Σ = Matrix(Hermitian(rand(100, 100) + 100I))
            μ = MvNormal(m, Σ)
            ν = MvNormal(m .+ τ, Σ)
            @test ot_cost(SqEuclidean(), μ, ν) ≈ norm(τ)^2

            x = rand(100, 10)
            T = ot_plan(SqEuclidean(), μ, ν)
            @test pdf(ν, mapslices(T, x; dims=1)) ≈ pdf(μ, x)
        end

        @testset "comparison to grid approximation" begin
            μ = MvNormal([0, 0], [1 0; 0 2])
            ν = MvNormal([10, 10], [2 0; 0 1])
            # Constructing circular grid approximation
            # Angular grid step
            θ = collect(0:0.2:(2π))
            θx = cos.(θ)
            θy = sin.(θ)
            # Radius grid step
            δ = collect(0:0.2:1)
            μsupp = [0.0 0.0]
            νsupp = [10.0 10.0]
            for i in δ[2:end]
                a = [θx .* i θy .* i * 2]
                b = [θx .* i * 2 θy .* i] .+ [10 10]
                μsupp = vcat(μsupp, a)
                νsupp = vcat(νsupp, b)
            end

            # Create discretized distribution
            μprobs = normalize!(pdf(μ, μsupp'), 1)
            νprobs = normalize!(pdf(ν, νsupp'), 1)
            C = pairwise(SqEuclidean(), μsupp', νsupp'; dims=2)
            for optimizer in (Tulip.Optimizer(), GLPK.Optimizer())
                @test emd2(μprobs, νprobs, C, optimizer) ≈ ot_cost(SqEuclidean(), μ, ν) rtol =
                    1e-3
            end

            # Use hcubature integration to perform ``\\int c(x,T(x)) d\\mu``
            T = ot_plan(SqEuclidean(), μ, ν)
            c_hcubature, _ = hcubature([-10, -10], [10, 10]) do x
                return sqeuclidean(x, T(x)) * pdf(μ, x)
            end
            @test ot_cost(SqEuclidean(), μ, ν) ≈ c_hcubature rtol = 1e-3
        end
    end
end
