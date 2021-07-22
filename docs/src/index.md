# MarkovBounds.jl 

[MarkovBounds.jl](https://github.com/FHoltorf/MarkovBounds.jl) is a Julia package seeking to automate the setup of moment bounding schemes for the analysis of jump-diffusion processes with the goal of enabling those unfamiliar with moment problems or optimization in general to apply moment bounding schemes regardless. To that end, MarkovBounds.jl automatically translates high-level problem data framing a jump-diffusion process as specified via convenient tools such as [Catalyst.jl](https://github.com/SciML/Catalyst.jl), [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)/[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) or [DynamicPolynomials.jl](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) into [sum-of-squares (SOS) programs](https://en.wikipedia.org/wiki/Sum-of-squares_optimization) and solve them via the existing optimization pipeline in Julia (see [SumOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl), [JuMP](https://github.com/jump-dev/JuMP.jl) and [MathOptInterface.jl](https://github.com/jump-dev/MathOptInterface.jl)). The solution of said SOS programs are returned to the user in form of theoretically guaranteed bounds on moments and other key statistics of the process under investigation.

![program structure](images\programstructure.PNG)

## When should you consider using moment bounding schemes?
Moment bounding schemes are limited by the capabilities of large-scale semidefinite programming. Given the current state-of-the-art, moment bounding schemes are practically limited to stochastic processes of low to medium dimensionality (< 10 states). Moreover, MarkovBounds.jl currently *only* supports processes in which the data can be fully characterized in terms of *polynomials* (see the [Background](@ref background) section for details). 

### Stationary Distributions
Moment bounding schemes have been found to perform remarkably well for the study of the stationary (or ergodic) statistics of stochastic processes. They have been found to provide high quality (often effectively tight) bounds on key statistics such as means, variances, Fano factors and more at a fraction of the computational time needed to provide similarly accurate sample statistics. 

### Error Control
The analysis of jump-diffusion processes relies traditionally on brute force sampling of paths of the process such that a sufficiently large sample of such paths exihibits similar statistics as the generating process. However, even though this approach works remarkably well over a wide range of applications and powerful software tools supporting this approach exist (see e.g. [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)), it may break down in certain settings. Most notably, the generation of sample paths can become prohibitively expensive when the process under investigation is highly stiff and/or volatile, rendering simulation expensive and/or convergence of sample statistics to the true statistics slow. In those cases, practitioners tend to fall back on approximation techniques such as ``\tau``-leaping, moment closure approximations and many more, which recover tractability at the cost of introducing unverifiable assumptions and an unknown error. Moment bounding schemes can be used to quantify/bound this error by providing hard bounds on key statistics, both in the transient and stationary setting. 

### Optimal Control Problems
Moment bounding schemes extend naturally to the application of (stochastic) optimal control problems with jump-diffusion processes embedded, where they allow to compute hard bounds on the optimal value. Such bounds can either be used in branch-and-bound algorithms for global optimization or simply as a way to certify optimality of a given control policy. Moreover, the bounding problems provide insights to control policy design by yielding a piecewise polynomial subsolution of the value function as byproduct. 

## Installation
You can install MarkovBounds.jl via Julia's package manager
```julia
 ]add https://github.com/FHoltorf/MarkovBounds.jl
```
Moreover, most functionalities of MarkovBounds.jl rely on a semidefinite programming (SDP) solver that is supported by JuMP/MathOptInterface. So please make sure that the SDP solver that you are planning on using is featured on this [list](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers). While you can choose any SDP solver from this [list](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers), we recommend to use either [Mosek](https://www.mosek.com/) or [SeDuMi](https://github.com/sqlp/sedumi) as they have shown the best results in our tests, both in terms of robustness and speed.

