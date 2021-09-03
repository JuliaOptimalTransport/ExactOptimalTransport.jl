# ExactOptimalTransport.jl <a href='https://juliaoptimaltransport.github.io/ExactOptimalTransport.jl/dev'><img src="docs/src/assets/logo.svg" align="right" height="138.5" /></a>

*Solving unregularized optimal transport problems with Julia*

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaOptimalTransport.github.io/ExactOptimalTransport.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaOptimalTransport.github.io/ExactOptimalTransport.jl/dev)
[![CI](https://github.com/JuliaOptimalTransport/ExactOptimalTransport.jl/workflows/CI/badge.svg?branch=main)](https://github.com/JuliaOptimalTransport/ExactOptimalTransport.jl/actions?query=workflow%3ACI+branch%3Amain)
[![Codecov](https://codecov.io/gh/JuliaOptimalTransport/ExactOptimalTransport.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaOptimalTransport/ExactOptimalTransport.jl)
[![Coveralls](https://coveralls.io/repos/github/JuliaOptimalTransport/ExactOptimalTransport.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaOptimalTransport/OptimalTransport.jl?branch=main)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

This package provides some [Julia](https://julialang.org/) implementations of algorithms for solving
unregularized optimal transport (Kantorovich) problems.

## Example

```julia
using ExactOptimalTransport
using Distances

# uniform histograms
μ = fill(1/250, 250)
ν = fill(1/200, 200)

# random cost matrix
C = pairwise(SqEuclidean(), rand(1, 250), rand(1, 200); dims=2)

# solve unregularized optimal transport problem
emd(μ, ν, C, ε)
```

Please see the documentation pages for further information.

## Related packages

- [OptimalTransport.jl](https://github.com/JuliaOptimalTransport/OptimalTransport.jl): Julia implementation of
algorithms for regularized optimal transport problems with GPU support.
- [StochasticOptimalTransport.jl](https://github.com/JuliaOptimalTransport/StochasticOptimalTransport.jl): Julia implementation of stochastic optimization algorithms for large-scale optimal transport.
- [PythonOT.jl](https://github.com/JuliaOptimalTransport/PythonOT.jl): Julia interface for the [Python Optimal Transport (POT) package](https://pythonot.github.io/).

## Contributing

Contributions are more than welcome! Please feel free to submit an issue or pull request in this repository.

## Note

This package was originally part of [OptimalTransport.jl](https://github.com/JuliaOptimalTransport/OptimalTransport.jl).


