using ExactOptimalTransport
using SafeTestsets

using Test

@testset "ExactOptimalTransport" begin
    @safetestset "Utilities" begin
        include("utils.jl")
    end

    @safetestset "Exact OT" begin
        include("exact.jl")
    end

    @safetestset "Wasserstein distance" begin
        include("wasserstein.jl")
    end

    @safetestset "Bures distance" begin
        include("bures.jl")
    end
end
