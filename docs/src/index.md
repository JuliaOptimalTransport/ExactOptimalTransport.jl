# ExactOptimalTransport.jl

ExactOptimalTransport.jl is a Julia package for solving the unregularized
optimal transport (Kantorovich) problems.

```@docs
emd
emd2
```

```@docs
ot_plan
ot_plan(::Any, ::ExactOptimalTransport.ContinuousUnivariateDistribution, ::ExactOptimalTransport.UnivariateDistribution)
ot_plan(::Any, ::ExactOptimalTransport.DiscreteNonParametric, ::ExactOptimalTransport.DiscreteNonParametric)
ot_plan(::ExactOptimalTransport.SqEuclidean, ::ExactOptimalTransport.Normal, ::ExactOptimalTransport.Normal)
ot_plan(::ExactOptimalTransport.SqEuclidean, ::ExactOptimalTransport.MvNormal, ::ExactOptimalTransport.MvNormal)
```

```@docs
ot_cost
ot_cost(::Any, ::ExactOptimalTransport.ContinuousUnivariateDistribution, ::ExactOptimalTransport.UnivariateDistribution)
ot_cost(::Any, ::ExactOptimalTransport.DiscreteNonParametric, ::ExactOptimalTransport.DiscreteNonParametric)
ot_cost(::ExactOptimalTransport.SqEuclidean, ::ExactOptimalTransport.Normal, ::ExactOptimalTransport.Normal)
ot_cost(::ExactOptimalTransport.SqEuclidean, ::ExactOptimalTransport.MvNormal, ::ExactOptimalTransport.MvNormal)
```

```@docs
wasserstein
squared2wasserstein
```
