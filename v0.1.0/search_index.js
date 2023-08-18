var documenterSearchIndex = {"docs":
[{"location":"#ExactOptimalTransport.jl","page":"Home","title":"ExactOptimalTransport.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ExactOptimalTransport.jl is a Julia package for solving the unregularized optimal transport (Kantorovich) problems.","category":"page"},{"location":"","page":"Home","title":"Home","text":"emd\nemd2","category":"page"},{"location":"#ExactOptimalTransport.emd","page":"Home","title":"ExactOptimalTransport.emd","text":"emd(μ, ν, C, optimizer)\n\nCompute the optimal transport plan γ for the Monge-Kantorovich problem with source histogram μ, target histogram ν, and cost matrix C of size (length(μ), length(ν)) which solves\n\ninf_γ  Π(μ ν) langle γ C rangle\n\nThe corresponding linear programming problem is solved with the user-provided optimizer. Possible choices are Tulip.Optimizer() and Clp.Optimizer() in the Tulip and Clp packages, respectively.\n\n\n\n\n\n","category":"function"},{"location":"#ExactOptimalTransport.emd2","page":"Home","title":"ExactOptimalTransport.emd2","text":"emd2(μ, ν, C, optimizer; plan=nothing)\n\nCompute the optimal transport cost (a scalar) for the Monge-Kantorovich problem with source histogram μ, target histogram ν, and cost matrix C of size (length(μ), length(ν)) which is given by\n\ninf_γ  Π(μ ν) langle γ C rangle\n\nThe corresponding linear programming problem is solved with the user-provided optimizer. Possible choices are Tulip.Optimizer() and Clp.Optimizer() in the Tulip and Clp packages, respectively.\n\nA pre-computed optimal transport plan may be provided.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"ot_plan\not_plan(::Any, ::ExactOptimalTransport.ContinuousUnivariateDistribution, ::ExactOptimalTransport.UnivariateDistribution)\not_plan(::Any, ::ExactOptimalTransport.DiscreteNonParametric, ::ExactOptimalTransport.DiscreteNonParametric)\not_plan(::ExactOptimalTransport.SqEuclidean, ::ExactOptimalTransport.Normal, ::ExactOptimalTransport.Normal)\not_plan(::ExactOptimalTransport.SqEuclidean, ::ExactOptimalTransport.MvNormal, ::ExactOptimalTransport.MvNormal)","category":"page"},{"location":"#ExactOptimalTransport.ot_plan","page":"Home","title":"ExactOptimalTransport.ot_plan","text":"ot_plan(c, μ, ν; kwargs...)\n\nCompute the optimal transport plan for the Monge-Kantorovich problem with source and target marginals μ and ν and cost c.\n\nThe optimal transport plan solves\n\ninf_gamma in Pi(mu nu) int c(x y)  mathrmdgamma(x y)\n\nwhere Pi(mu nu) denotes the couplings of mu and nu.\n\nSee also: ot_cost\n\n\n\n\n\n","category":"function"},{"location":"#ExactOptimalTransport.ot_plan-Tuple{Any, Distributions.Distribution{Distributions.Univariate, Distributions.Continuous}, Distributions.UnivariateDistribution{S} where S<:Distributions.ValueSupport}","page":"Home","title":"ExactOptimalTransport.ot_plan","text":"ot_plan(c, μ::ContinuousUnivariateDistribution, ν::UnivariateDistribution)\n\nCompute the optimal transport plan for the Monge-Kantorovich problem with univariate distributions μ and ν as source and target marginals and cost function c of the form c(x y) = h(x - y) where h is a convex function.\n\nIn this setting, the optimal transport plan is the Monge map\n\nT = F_nu^-1 circ F_mu\n\nwhere F_mu is the cumulative distribution function of μ and F_nu^-1 is the quantile function of ν.\n\nSee also: ot_cost, emd\n\n\n\n\n\n","category":"method"},{"location":"#ExactOptimalTransport.ot_plan-Tuple{Any, Distributions.DiscreteNonParametric, Distributions.DiscreteNonParametric}","page":"Home","title":"ExactOptimalTransport.ot_plan","text":"ot_plan(c, μ::DiscreteNonParametric, ν::DiscreteNonParametric)\n\nCompute the optimal transport cost for the Monge-Kantorovich problem with univariate discrete distributions μ and ν as source and target marginals and cost function c of the form c(x y) = h(x - y) where h is a convex function.\n\nIn this setting, the optimal transport plan can be computed analytically. It is returned as a sparse matrix.\n\nSee also: ot_cost, emd\n\n\n\n\n\n","category":"method"},{"location":"#ExactOptimalTransport.ot_plan-Tuple{Distances.SqEuclidean, Distributions.Normal, Distributions.Normal}","page":"Home","title":"ExactOptimalTransport.ot_plan","text":"ot_plan(::SqEuclidean, μ::Normal, ν::Normal)\n\nCompute the optimal transport plan for the Monge-Kantorovich problem with normal distributions μ and ν as source and target marginals and cost function c(x y) = x - y_2^2.\n\nSee also: ot_cost, emd\n\n\n\n\n\n","category":"method"},{"location":"#ExactOptimalTransport.ot_plan-Tuple{Distances.SqEuclidean, Distributions.MvNormal, Distributions.MvNormal}","page":"Home","title":"ExactOptimalTransport.ot_plan","text":"ot_plan(::SqEuclidean, μ::MvNormal, ν::MvNormal)\n\nCompute the optimal transport plan for the Monge-Kantorovich problem with multivariate normal distributions μ and ν as source and target marginals and cost function c(x y) = x - y_2^2.\n\nIn this setting, for mu = mathcalN(m_mu Sigma_mu) and nu = mathcalN(m_nu Sigma_nu), the optimal transport plan is the Monge map\n\nT colon x mapsto m_nu\n+ Sigma_mu^-12\nbig(Sigma_mu^12 Sigma_nu Sigma_mu^12big)^12Sigma_mu^-12\n(x - m_mu)\n\nSee also: ot_cost, emd\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"ot_cost\not_cost(::Any, ::ExactOptimalTransport.ContinuousUnivariateDistribution, ::ExactOptimalTransport.UnivariateDistribution)\not_cost(::Any, ::ExactOptimalTransport.DiscreteNonParametric, ::ExactOptimalTransport.DiscreteNonParametric)\not_cost(::ExactOptimalTransport.SqEuclidean, ::ExactOptimalTransport.Normal, ::ExactOptimalTransport.Normal)\not_cost(::ExactOptimalTransport.SqEuclidean, ::ExactOptimalTransport.MvNormal, ::ExactOptimalTransport.MvNormal)","category":"page"},{"location":"#ExactOptimalTransport.ot_cost","page":"Home","title":"ExactOptimalTransport.ot_cost","text":"ot_cost(c, μ, ν; kwargs...)\n\nCompute the optimal transport cost for the Monge-Kantorovich problem with source and target marginals μ and ν and cost c.\n\nThe optimal transport cost is the scalar value\n\ninf_gamma in Pi(mu nu) int c(x y)  mathrmdgamma(x y)\n\nwhere Pi(mu nu) denotes the couplings of mu and nu.\n\nSee also: ot_plan\n\n\n\n\n\n","category":"function"},{"location":"#ExactOptimalTransport.ot_cost-Tuple{Any, Distributions.Distribution{Distributions.Univariate, Distributions.Continuous}, Distributions.UnivariateDistribution{S} where S<:Distributions.ValueSupport}","page":"Home","title":"ExactOptimalTransport.ot_cost","text":"ot_cost(\n    c, μ::ContinuousUnivariateDistribution, ν::UnivariateDistribution; plan=nothing\n)\n\nCompute the optimal transport cost for the Monge-Kantorovich problem with univariate distributions μ and ν as source and target marginals and cost function c of the form c(x y) = h(x - y) where h is a convex function.\n\nIn this setting, the optimal transport cost can be computed as\n\nint_0^1 c(F_mu^-1(x) F_nu^-1(x)) mathrmdx\n\nwhere F_mu^-1 and F_nu^-1 are the quantile functions of μ and ν, respectively.\n\nA pre-computed optimal transport plan may be provided.\n\nSee also: ot_plan, emd2\n\n\n\n\n\n","category":"method"},{"location":"#ExactOptimalTransport.ot_cost-Tuple{Any, Distributions.DiscreteNonParametric, Distributions.DiscreteNonParametric}","page":"Home","title":"ExactOptimalTransport.ot_cost","text":"ot_cost(\n    c, μ::DiscreteNonParametric, ν::DiscreteNonParametric; plan=nothing\n)\n\nCompute the optimal transport cost for the Monge-Kantorovich problem with discrete univariate distributions μ and ν as source and target marginals and cost function c of the form c(x y) = h(x - y) where h is a convex function.\n\nIn this setting, the optimal transport cost can be computed analytically.\n\nA pre-computed optimal transport plan may be provided.\n\nSee also: ot_plan, emd2\n\n\n\n\n\n","category":"method"},{"location":"#ExactOptimalTransport.ot_cost-Tuple{Distances.SqEuclidean, Distributions.Normal, Distributions.Normal}","page":"Home","title":"ExactOptimalTransport.ot_cost","text":"ot_cost(::SqEuclidean, μ::Normal, ν::Normal)\n\nCompute the squared 2-Wasserstein distance between univariate normal distributions μ and ν as source and target marginals.\n\nSee also: ot_plan, emd2\n\n\n\n\n\n","category":"method"},{"location":"#ExactOptimalTransport.ot_cost-Tuple{Distances.SqEuclidean, Distributions.MvNormal, Distributions.MvNormal}","page":"Home","title":"ExactOptimalTransport.ot_cost","text":"ot_cost(::SqEuclidean, μ::MvNormal, ν::MvNormal)\n\nCompute the squared 2-Wasserstein distance between normal distributions μ and ν as source and target marginals.\n\nIn this setting, the optimal transport cost can be computed as\n\nW_2^2(mu nu) = m_mu - m_nu ^2 + mathcalB(Sigma_mu Sigma_nu)^2\n\nwhere mu = mathcalN(m_mu Sigma_mu), nu = mathcalN(m_nu Sigma_nu), and mathcalB is the Bures metric.\n\nSee also: ot_plan, emd2\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"wasserstein\nsquared2wasserstein","category":"page"},{"location":"#ExactOptimalTransport.wasserstein","page":"Home","title":"ExactOptimalTransport.wasserstein","text":"wasserstein(μ, ν; metric=Euclidean(), p=Val(1), kwargs...)\n\nCompute the p-Wasserstein distance with respect to the metric between measures μ and ν.\n\nOrder p can be provided as a scalar of type Real or as a parameter of a value type Val(p). For certain combinations of metric and p, such as metric=Euclidean() and p=Val(2), the computations are more efficient if p is specified as a value type. The remaining keyword arguments are forwarded to ot_cost.\n\nSee also: squared2wasserstein, ot_cost\n\n\n\n\n\n","category":"function"},{"location":"#ExactOptimalTransport.squared2wasserstein","page":"Home","title":"ExactOptimalTransport.squared2wasserstein","text":"squared2wasserstein(μ, ν; metric=Euclidean(), kwargs...)\n\nCompute the squared 2-Wasserstein distance with respect to the metric between measures μ and ν.\n\nThe remaining keyword arguments are forwarded to ot_cost.\n\nSee also: wasserstein, ot_cost\n\n\n\n\n\n","category":"function"}]
}
