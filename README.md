# MarkovBounds

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://FHoltorf.github.io/MarkovBounds.jl/dev)
[![Build Status](https://github.com/FHoltorf/MarkovBounds.jl/workflows/CI/badge.svg)](https://github.com/FHoltorf/MarkovBounds.jl/actions)
[![Coverage](https://codecov.io/gh/FHoltorf/MarkovBounds.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/FHoltorf/MarkovBounds.jl)

MarkovBounds.jl is a meta-package based on SumOfSquares.jl which is intended to simplify the use of moment bounding schemes for the analysis of jump-diffusion processes with polynomial data and associated optimal control problems. The hope is that this package will enable non-expert users, with special emphasis on the stochastic chemical kinetics community, to apply these techniques by automating the translation of high-level problem data specified by the user through symbolic tools (Catalyst.jl, Symbolics.jl or DynamicPolynomials.jl) into moment bounding problems. The problems can the be solved using the existing optimization infrastructure provided by the MathOptInterface. 

Moment bounding schemes enable the computation of *hard, theoretically guaranteed* bounds on moments and many related statistics associated with jump-diffusion processes, both for the transient case and the stationary distribution of the process. These statistics include for example
* means
* variances
* the volume of confidence ellipsoids of polynomial expressions
* the probability mass of key events characterized by membership of the process state in a basic closed semialgebraic set
These techniques further extend naturally to computing bounds on the optimal value of many finite horizon and (discounted) infinite horizon optimal control problems. 
