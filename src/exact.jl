"""
    ot_plan(c, μ, ν; kwargs...)

Compute the optimal transport plan for the Monge-Kantorovich problem with source and target
marginals `μ` and `ν` and cost `c`.

The optimal transport plan solves
```math
\\inf_{\\gamma \\in \\Pi(\\mu, \\nu)} \\int c(x, y) \\, \\mathrm{d}\\gamma(x, y)
```
where ``\\Pi(\\mu, \\nu)`` denotes the couplings of ``\\mu`` and ``\\nu``.

See also: [`ot_cost`](@ref)
"""
function ot_plan end

"""
    ot_cost(c, μ, ν; kwargs...)

Compute the optimal transport cost for the Monge-Kantorovich problem with source and target
marginals `μ` and `ν` and cost `c`.

The optimal transport cost is the scalar value
```math
\\inf_{\\gamma \\in \\Pi(\\mu, \\nu)} \\int c(x, y) \\, \\mathrm{d}\\gamma(x, y)
```
where ``\\Pi(\\mu, \\nu)`` denotes the couplings of ``\\mu`` and ``\\nu``.

See also: [`ot_plan`](@ref)
"""
function ot_cost end

#############
# Discrete OT
#############

"""
    emd(μ, ν, C, optimizer)

Compute the optimal transport plan `γ` for the Monge-Kantorovich problem with source
histogram `μ`, target histogram `ν`, and cost matrix `C` of size `(length(μ), length(ν))`
which solves
```math
\\inf_{γ ∈ Π(μ, ν)} \\langle γ, C \\rangle.
```

The corresponding linear programming problem is solved with the user-provided `optimizer`.
Possible choices are `Tulip.Optimizer()` and `Clp.Optimizer()` in the `Tulip` and `Clp`
packages, respectively.
"""
function emd(μ, ν, C, model::MOI.ModelLike)
    # check size of cost matrix
    nμ = length(μ)
    nν = length(ν)
    size(C) == (nμ, nν) || error("cost matrix `C` must be of size `(length(μ), length(ν))`")
    nC = length(C)

    # define variables
    x = MOI.add_variables(model, nC)
    xmat = reshape(x, nμ, nν)

    # define objective function
    T = float(eltype(C))
    zero_T = zero(T)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(float.(vec(C)), x), zero_T),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    # add non-negativity constraints
    for xi in x
        MOI.add_constraint(model, xi, MOI.GreaterThan(zero_T))
    end

    # add constraints for source
    for (xrow, μi) in zip(eachrow(xmat), μ)
        f = MOI.ScalarAffineFunction(
            [MOI.ScalarAffineTerm(one(μi), xi) for xi in xrow], zero(μi)
        )
        MOI.add_constraint(model, f, MOI.EqualTo(μi))
    end

    # add constraints for target
    for (xcol, νi) in zip(eachcol(xmat), ν)
        f = MOI.ScalarAffineFunction(
            [MOI.ScalarAffineTerm(one(νi), xi) for xi in xcol], zero(νi)
        )
        MOI.add_constraint(model, f, MOI.EqualTo(νi))
    end

    # compute optimal solution
    MOI.optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    status === MOI.OPTIMAL || error("failed to compute optimal transport plan: ", status)
    p = MOI.get(model, MOI.VariablePrimal(), x)
    γ = reshape(p, nμ, nν)

    return γ
end

"""
    emd2(μ, ν, C, optimizer; plan=nothing)

Compute the optimal transport cost (a scalar) for the Monge-Kantorovich problem with source
histogram `μ`, target histogram `ν`, and cost matrix `C` of size `(length(μ), length(ν))`
which is given by
```math
\\inf_{γ ∈ Π(μ, ν)} \\langle γ, C \\rangle.
```

The corresponding linear programming problem is solved with the user-provided `optimizer`.
Possible choices are `Tulip.Optimizer()` and `Clp.Optimizer()` in the `Tulip` and `Clp`
packages, respectively.

A pre-computed optimal transport `plan` may be provided.
"""
function emd2(μ, ν, C, optimizer; plan=nothing)
    γ = if plan === nothing
        # compute optimal transport plan
        emd(μ, ν, C, optimizer)
    else
        # check dimensions
        size(C) == (length(μ), length(ν)) ||
            error("cost matrix `C` must be of size `(length(μ), length(ν))`")
        size(plan) == size(C) || error(
            "optimal transport plan `plan` and cost matrix `C` must be of the same size",
        )
        plan
    end
    return dot(γ, C)
end

###################################
# Semidiscrete and continuous 1D OT
###################################

"""
    ot_plan(c, μ::ContinuousUnivariateDistribution, ν::UnivariateDistribution)

Compute the optimal transport plan for the Monge-Kantorovich problem with univariate
distributions `μ` and `ν` as source and target marginals and cost function `c` of
the form ``c(x, y) = h(|x - y|)`` where ``h`` is a convex function.

In this setting, the optimal transport plan is the Monge map
```math
T = F_\\nu^{-1} \\circ F_\\mu
```
where ``F_\\mu`` is the cumulative distribution function of `μ` and ``F_\\nu^{-1}`` is the
quantile function of `ν`.

See also: [`ot_cost`](@ref), [`emd`](@ref)
"""
function ot_plan(c, μ::ContinuousUnivariateDistribution, ν::UnivariateDistribution)
    # Use T instead of γ to indicate that this is a Monge map.
    T(x) = quantile(ν, cdf(μ, x))
    return T
end

"""
    ot_cost(
        c, μ::ContinuousUnivariateDistribution, ν::UnivariateDistribution; plan=nothing
    )

Compute the optimal transport cost for the Monge-Kantorovich problem with univariate
distributions `μ` and `ν` as source and target marginals and cost function `c` of
the form ``c(x, y) = h(|x - y|)`` where ``h`` is a convex function.

In this setting, the optimal transport cost can be computed as
```math
\\int_0^1 c(F_\\mu^{-1}(x), F_\\nu^{-1}(x)) \\mathrm{d}x
```
where ``F_\\mu^{-1}`` and ``F_\\nu^{-1}`` are the quantile functions of `μ` and `ν`,
respectively.

A pre-computed optimal transport `plan` may be provided.

See also: [`ot_plan`](@ref), [`emd2`](@ref)
"""
function ot_cost(
    c, μ::ContinuousUnivariateDistribution, ν::UnivariateDistribution; plan=nothing
)
    cost, _ = if plan === nothing
        quadgk(0, 1) do q
            return c(quantile(μ, q), quantile(ν, q))
        end
    else
        quadgk(0, 1) do q
            x = quantile(μ, q)
            return c(x, plan(x))
        end
    end
    return cost
end

################
# Discrete 1D OT
################

# internal iterator for discrete one-dimensional OT problems
# it returns tuples that consist of the indices of the source and target histograms
# and the optimal flow between the corresponding points
struct Discrete1DOTIterator{T,M,N}
    mu::M
    nu::N
end

# histograms `μ` and `ν` are expected to be iterators of the histograms where the
# corresponding support is sorted
function Discrete1DOTIterator(μ, ν)
    T = Base.promote_eltype(μ, ν)
    return Discrete1DOTIterator{T,typeof(μ),typeof(ν)}(μ, ν)
end

Base.IteratorEltype(::Type{<:Discrete1DOTIterator}) = Base.HasEltype()
Base.eltype(::Type{<:Discrete1DOTIterator{T}}) where {T} = Tuple{Int,Int,T}

Base.length(d::Discrete1DOTIterator) = length(d.mu) + length(d.nu) - 1

# we iterate through the source and target histograms
function Base.iterate(
    d::Discrete1DOTIterator{T}, (i, j, μnext, νnext)=(1, 1, iterate(d.mu), iterate(d.nu))
) where {T}
    # if we are done with iterating through the source and/or target histogram,
    # iteration is stopped
    if μnext === nothing || νnext === nothing
        return nothing
    end

    # unpack next values and states of the source and target histograms
    μiter, μstate = μnext
    νiter, νstate = νnext

    # compute next value of the iterator: indices of source and target histograms
    # and optimal flow between the corresponding points
    min_iter, max_iter = minmax(μiter, νiter)
    iter = (i, j, min_iter)

    # compute next state of the iterator
    diff = max_iter - min_iter
    state = if μiter < max_iter
        # move forward in the source histogram
        (i + 1, j, iterate(d.mu, μstate), (diff, νstate))
    else
        # move forward in the target histogram
        (i, j + 1, (diff, μstate), iterate(d.nu, νstate))
    end

    return iter, state
end

"""
    ot_plan(c, μ::DiscreteNonParametric, ν::DiscreteNonParametric)

Compute the optimal transport cost for the Monge-Kantorovich problem with univariate
discrete distributions `μ` and `ν` as source and target marginals and cost function `c`
of the form ``c(x, y) = h(|x - y|)`` where ``h`` is a convex function.

In this setting, the optimal transport plan can be computed analytically. It is returned as
a sparse matrix.

See also: [`ot_cost`](@ref), [`emd`](@ref)
"""
function ot_plan(_, μ::DiscreteNonParametric, ν::DiscreteNonParametric)
    # Unpack the probabilities of the two distributions
    # Note: support of `DiscreteNonParametric` is sorted
    μprobs = probs(μ)
    νprobs = probs(ν)
    T = Base.promote_eltype(μprobs, νprobs)

    return if μprobs isa FillArrays.AbstractFill &&
        νprobs isa FillArrays.AbstractFill &&
        length(μprobs) == length(νprobs)
        # Special case: discrete uniform distributions of the same "size"
        k = length(μprobs)
        sparse(1:k, 1:k, T(first(μprobs)), k, k)
    else
        # Generic case
        # Create the iterator
        iter = Discrete1DOTIterator(μprobs, νprobs)

        # create arrays for the indices of the two histograms and the optimal flow between the
        # corresponding points
        n = length(iter)
        I = Vector{Int}(undef, n)
        J = Vector{Int}(undef, n)
        W = Vector{T}(undef, n)

        # compute the sparse optimal transport plan
        @inbounds for (idx, (i, j, w)) in enumerate(iter)
            I[idx] = i
            J[idx] = j
            W[idx] = w
        end
        sparse(I, J, W, length(μprobs), length(νprobs))
    end
end

"""
    ot_cost(
        c, μ::DiscreteNonParametric, ν::DiscreteNonParametric; plan=nothing
    )

Compute the optimal transport cost for the Monge-Kantorovich problem with discrete
univariate distributions `μ` and `ν` as source and target marginals and cost function `c`
of the form ``c(x, y) = h(|x - y|)`` where ``h`` is a convex function.

In this setting, the optimal transport cost can be computed analytically.

A pre-computed optimal transport `plan` may be provided.

See also: [`ot_plan`](@ref), [`emd2`](@ref)
"""
function ot_cost(c, μ::DiscreteNonParametric, ν::DiscreteNonParametric; plan=nothing)
    # Extract support and probabilities of discrete distributions
    # Note: support of `DiscreteNonParametric` is sorted
    μsupport = support(μ)
    νsupport = support(ν)
    μprobs = probs(μ)
    νprobs = probs(ν)

    return if μprobs isa FillArrays.AbstractFill &&
        νprobs isa FillArrays.AbstractFill &&
        length(μprobs) == length(νprobs)
        # Special case: discrete uniform distributions of the same "size"
        # In this case we always just compute `sum(c.(μsupport .- νsupport))` and scale it
        # We use pairwise summation and avoid allocations
        # (https://github.com/JuliaLang/julia/pull/31020)
        T = Base.promote_eltype(μprobs, νprobs)
        T(first(μprobs)) *
        sum(Broadcast.instantiate(Broadcast.broadcasted(c, μsupport, νsupport)))
    else
        # Generic case 
        _ot_cost(c, μsupport, μprobs, νsupport, νprobs, plan)
    end
end

# compute cost from scratch if no plan is provided
function _ot_cost(c, μsupport, μprobs, νsupport, νprobs, ::Nothing)
    # create the iterator
    iter = Discrete1DOTIterator(μprobs, νprobs)

    # compute the cost
    return sum(w * c(μsupport[i], νsupport[j]) for (i, j, w) in iter)
end

# if a sparse plan is provided, we just iterate through the non-zero entries
function _ot_cost(c, μsupport, _, νsupport, _, plan::SparseMatrixCSC)
    # extract non-zero flows
    I, J, W = findnz(plan)

    # compute the cost
    return sum(w * c(μsupport[i], νsupport[j]) for (i, j, w) in zip(I, J, W))
end

# fallback: compute cost matrix (probably often faster to compute cost from scratch)
function _ot_cost(c, μsupport, _, νsupport, _, plan)
    return dot(plan, StatsBase.pairwise(c, μsupport, νsupport))
end

################
# OT Gaussians
################

"""
    ot_cost(::SqEuclidean, μ::MvNormal, ν::MvNormal)

Compute the squared 2-Wasserstein distance between normal distributions `μ` and `ν` as
source and target marginals.

In this setting, the optimal transport cost can be computed as
```math
W_2^2(\\mu, \\nu) = \\|m_\\mu - m_\\nu \\|^2 + \\mathcal{B}(\\Sigma_\\mu, \\Sigma_\\nu)^2,
```
where ``\\mu = \\mathcal{N}(m_\\mu, \\Sigma_\\mu)``,
``\\nu = \\mathcal{N}(m_\\nu, \\Sigma_\\nu)``, and ``\\mathcal{B}`` is the Bures metric.

See also: [`ot_plan`](@ref), [`emd2`](@ref)
"""
function ot_cost(::SqEuclidean, μ::MvNormal, ν::MvNormal)
    return sqeuclidean(μ.μ, ν.μ) + sqbures(μ.Σ, ν.Σ)
end

"""
    ot_cost(::SqEuclidean, μ::Normal, ν::Normal)

Compute the squared 2-Wasserstein distance between univariate normal distributions `μ` and
`ν` as source and target marginals.

See also: [`ot_plan`](@ref), [`emd2`](@ref)
"""
function ot_cost(::SqEuclidean, μ::Normal, ν::Normal)
    return (μ.μ - ν.μ)^2 + (μ.σ - ν.σ)^2
end

"""
    ot_plan(::SqEuclidean, μ::MvNormal, ν::MvNormal)

Compute the optimal transport plan for the Monge-Kantorovich problem with multivariate
normal distributions `μ` and `ν` as source and target marginals and cost function
``c(x, y) = \\|x - y\\|_2^2``.

In this setting, for ``\\mu = \\mathcal{N}(m_\\mu, \\Sigma_\\mu)`` and
``\\nu = \\mathcal{N}(m_\\nu, \\Sigma_\\nu)``, the optimal transport plan is the Monge
map
```math
T \\colon x \\mapsto m_\\nu
+ \\Sigma_\\mu^{-1/2}
{\\big(\\Sigma_\\mu^{1/2} \\Sigma_\\nu \\Sigma_\\mu^{1/2}\\big)}^{1/2}\\Sigma_\\mu^{-1/2}
(x - m_\\mu).
```

See also: [`ot_cost`](@ref), [`emd`](@ref)
"""
function ot_plan(::SqEuclidean, μ::MvNormal, ν::MvNormal)
    Σμsqrt = μ.Σ^(-1 / 2)
    A = Σμsqrt * sqrt(_gaussian_ot_A(μ.Σ, ν.Σ)) * Σμsqrt
    mμ = μ.μ
    mν = ν.μ
    T(x) = mν + A * (x - mμ)
    return T
end

"""
    ot_plan(::SqEuclidean, μ::Normal, ν::Normal)

Compute the optimal transport plan for the Monge-Kantorovich problem with
normal distributions `μ` and `ν` as source and target marginals and cost function
``c(x, y) = \\|x - y\\|_2^2``.

See also: [`ot_cost`](@ref), [`emd`](@ref)
"""
function ot_plan(::SqEuclidean, μ::Normal, ν::Normal)
    mμ = μ.μ
    mν = ν.μ
    a = ν.σ / μ.σ
    T(x) = mν + a * (x - mμ)
    return T
end
