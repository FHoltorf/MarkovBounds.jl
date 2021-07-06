# MarkovBounds

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://FHoltorf.github.io/MarkovBounds.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://FHoltorf.github.io/MarkovBounds.jl/dev)
[![Build Status](https://github.com/FHoltorf/MarkovBounds.jl/workflows/CI/badge.svg)](https://github.com/FHoltorf/MarkovBounds.jl/actions)
[![Coverage](https://codecov.io/gh/FHoltorf/MarkovBounds.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/FHoltorf/MarkovBounds.jl)

MarkovBounds.jl is a meta-package based on SumOfSquares.jl which is intended to simplify the use of moment bounding schemes for the analysis of jump-diffusion processes with polynomial data and associated optimal control problems via convex optimization. The hope is that this package will enable non-expert users, with special emphasis on the stochastic chemical kinetics community, to apply these techniques by automating problem setup and formulation such that the user only is required to provide high-level problem data such as a reaction network for example. 

Moment bounding schemes enable the computation of *hard, theoretically guaranteed* bounds on moments and many related statistics such as variances, Fano factors, or confidence ellipsoid volumes associated with jump-diffusion processes. Likewise, bounds on the optimal value of many finite and (discounted) infinite horizon optimal control problems can be computed efficiently.
