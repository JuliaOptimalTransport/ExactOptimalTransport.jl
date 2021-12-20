module ExactOptimalTransport

using Distances
using MathOptInterface
using Distributions
using FillArrays
using PDMats
using QuadGK
using StatsBase: StatsBase

using LinearAlgebra
using SparseArrays

export emd, emd2
export ot_cost, ot_plan, wasserstein, squared2wasserstein
export discretemeasure

const MOI = MathOptInterface

include("distances/bures.jl")
include("utils.jl")
include("exact.jl")
include("wasserstein.jl")

end
