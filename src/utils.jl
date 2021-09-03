struct FiniteDiscreteMeasure{X<:AbstractVector,P<:AbstractVector}
    support::X
    p::P

    function FiniteDiscreteMeasure{X,P}(support::X, p::P) where {X,P}
        length(support) == length(p) || error("length of `support` and `p` must be equal")
        isprobvec(p) || error("`p` must be a probability vector")
        return new{X,P}(support, p)
    end
end

"""
    discretemeasure(
        support::AbstractVector,
        probs::AbstractVector{<:Real}=fill(inv(length(support)), length(support))
    )

Construct a finite discrete probability measure with `support` and corresponding
`probabilities`. If the probability vector argument is not passed, then
equal probability is assigned to each entry in the support.

# Examples
```julia
using KernelFunctions
# rows correspond to samples
μ = discretemeasure(RowVecs(rand(7,3)), normalize!(rand(10),1))

# columns correspond to samples, each with equal probability
ν = discretemeasure(ColVecs(rand(3,12)))
```

!!! note
    If `support` is a 1D vector, the constructed measure will be sorted,
    e.g. for `mu = discretemeasure([3, 1, 2],[0.5, 0.2, 0.3])`, then
    `mu.support` will be `[1, 2, 3]` and `mu.p` will be `[0.2, 0.3, 0.5]`.
"""
function discretemeasure(
    support::AbstractVector{<:Real},
    probs::AbstractVector{<:Real}=fill(inv(length(support)), length(support)),
)
    return DiscreteNonParametric(support, probs)
end
function discretemeasure(
    support::AbstractVector,
    probs::AbstractVector{<:Real}=fill(inv(length(support)), length(support)),
)
    return FiniteDiscreteMeasure{typeof(support),typeof(probs)}(support, probs)
end

Distributions.support(d::FiniteDiscreteMeasure) = d.support
Distributions.probs(d::FiniteDiscreteMeasure) = d.p
