```@meta
CurrentModule = MarkovBounds
```

# MarkovBounds.jl 

[MarkovBounds.jl](https://github.com/FHoltorf/MarkovBounds.jl) is a Julia package seeking to automate the setup of moment bounding schemes for the analysis of jump-diffusion processes with the goal of enabling those unfamiliar with moment problems or optimization in general to apply moment bounding schemes regardless. To that end, MarkovBounds.jl automatically translates high-level problem data framing a jump-diffusion process as specified via convenient tools such as [Catalyst.jl](https://github.com/SciML/Catalyst.jl), [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)/[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) or [DynamicPolynomials.jl](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) into [sum-of-squares (SOS) programs](https://en.wikipedia.org/wiki/Sum-of-squares_optimization) and solve them via the existing optimization pipeline in Julia (see [SumOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl), [JuMP](https://github.com/jump-dev/JuMP.jl) and [MathOptInterface.jl](https://github.com/jump-dev/MathOptInterface.jl)). The solution of said SOS programs are returned to the user in form of theoretically guaranteed bounds on moments and other key statistics of the process under investigation.

![program structure](C:\Users\chemegrad2018\.julia\dev\MarkovBounds\docs\src\images\programstructure.PNG)



# When should you consider using moment bounding schemes?
Moment bounding schemes are limited by the capabilities of large-scale semidefinite programming. Given the current state-of-the-art, moment bounding schemes are practically limited to stochastic processes of low to medium dimensionality (< 10 states). Moreover, MarkovBounds.jl currently *only* supports processes in which the data can be fully characterized in terms of *polynomials* (see [Background](#background-on-moment-bounding-schemes) for details). 

## Stationary Distributions
Moment bounding schemes have been found to perform remarkably well for the study of the stationary (or ergodic) statistics of stochastic processes. They have been found to provide high quality (often effectively tight) bounds on key statistics such as means, variances, Fano factors and more at a fraction of the computational time needed to provide similarly accurate sample statistics. 

## Error Control
The analysis of jump-diffusion processes relies traditionally on brute force sampling of paths of the process such that a sufficiently large sample of such paths exihibits similar statistics as the generating process. However, even though this approach works remarkably well over a wide range of applications and powerful software tools supporting this approach exist (see e.g. [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)), it may break down in certain settings. Most notably, the generation of sample paths can become prohibitively expensive when the process under investigation is highly stiff and/or volatile, rendering simulation expensive and/or convergence of sample statistics to the true statistics slow. In those cases, practitioners tend to fall back on approximation techniques such as $\tau$-leaping, moment closure approximations and many more, which recover tractability at the cost of introducing unverifiable assumptions and an unknown error. Moment bounding schemes can be used to quantify/bound this error by providing hard bounds on key statistics, both in the transient and stationary setting. 

## Optimal Control Problems
Moment bounding schemes extend naturally to the application of (stochastic) optimal control problems with jump-diffusion processes embedded, where they allow to compute hard bounds on the optimal value. Such bounds can either be used in branch-and-bound algorithms for global optimization or simply as a way to certify optimality of a given control policy. Moreover, the bounding problems provide insights to control policy design by yielding a piecewise polynomial subsolution of the value function as byproduct. 

# Background on Moment Bounding Schemes

## Jump-Diffusion Processes
A jump-diffusion process is dynamical system combining a deterministic evolution of the system state, called drift, with a stochastic component modeling stochastic vibrations, called diffusion, and another stochastic component modeling discrete changes, called jumps. The evolution of the process state $ x_t $ over time $ t $ through its state space $ X \subset \mathbb{R}^n $ is governed by the following stochastic differential equation
$$
    dx_t = f(x_t) \, dt + g(x_t) \, dW_t + \sum_{i=1}^{n_R} h(x_t) \, dN_{a_i(x_t)}
$$
where $W_t$ denotes a standard $\mathbb{R}^m$-Brownian motion and $N_{a_i(x_t)}$ a standard Poisson counter with rate $a_i(x)$. The problem data is considered

* drift coefficient $f:\mathbb{R}^n \to \mathbb{R}^n$
* diffusion matrix $gg^\top :  \mathbb{R}^n \to \mathbb{R}^{n\times n}$ (or diffusion coefficient $ g:\mathbb{R}^n \to \mathbb{R}^{n \times m} $)
* arrival rates $a_i : \mathbb{R}^n \to \mathbb{R}$, $i  = 1,\dots, n_R$
* jumps $h_i:\mathbb{R}^n \to \mathbb{R}^n$, $i  = 1,\dots, n_R$
* state space $X$ 

A fundamental assumption in MarkovBounds.jl (and to a large extent moment bounding schemes inherently) is that the data of the jump-diffusion process under investigation can be fully characterized in terms of polynomials, i.e., all functions listed above are polynomials (component-wise) and the state space is (or can at least be outer approximated by) a [basic closed semialgebraic set](https://www.mit.edu/~parrilo/cdc03_workshop/10_positivstellensatz_2003_12_07_02_screen.pdf). Throughout, we will refer to processes that satisfies this assumption as polynomial jump-diffusion processes. A wide range of problems, in particular in the realm of stochastic chemical kinetics, lend themselves to be modeled in terms of polynomial jump-diffusion processes; if this assumption, however, is not satisfied, a simple but limited workaround is to find approximations of the data in terms of polynomials and apply the moment bounding scheme in a second step.

## Moment Bounding Schemes
The core idea behind moment bounding schemes is rather simple. But to explain it, we first need to establish some notation: Let $y_i(t) = \mathbb{E} \left[ \prod_{k=1}^n x_k(t)^{i_k} \right]$ denote the moment corresponding to the multi-index $ i \in \mathbb{N}_{0}^n $ of a polynomial jump-diffusion process as defined the previous section. Similarly, let $ \mathbf{y} _ {q}(t) $  be the truncated sequence of all multivariate moments of the process up to order $q \in \mathbb{N}$, i.e., $\mathbf{y} _ q(t) = \{ y_i(t) | |i| \leq q \} $. Due to the notorious moment closure problem, $ \mathbf{y} _ q(t) $ cannot in general be computed directly via simple simulation. To circumvent this issue, moment bounding schemes seek to identify a proxy for $\mathbf{y}_q(t)$, say $ \tilde{\mathbf{y}}_q (t) $, which minimizes (or maximizes if upper bound is sought) a certain statistic of the process under investigation, while ensuring that $\tilde{\mathbf{y}}_q(t)$ remains in certain ways consistent with the process under investigation (we will see shortly what that means concretely). Slightly more formally, we seek to solve an optimization problem of the form
$$
\begin{align} \inf_{\tilde{\mathbf{y}}_q} \quad &\int_{0}^T l^\top \tilde{\mathbf{y}}_q(t) \, dt + m^\top \tilde{\mathbf{y}}_q(T) \\
\text{s.t.} \quad & \tilde{\mathbf{y}}_q \text{ satisfies necessary consistency conditions.} \end{align}
$$
The key insight underpinning all moment bounding schemes now is that a suitable choice of "necessary consistency conditions" turns the above "pseudo" optimization problem into a convex optimization problem known as generalized moment problem. The practical value of this observations lies in the fact that strong convex relaxations of these generalized moment problems are easily constructed and they can be readily solved with off-the-shelve semidefinite programming (SDP) solvers such as Mosek, SeDuMi or SDPT3. 

But what are these "necessary consistency conditions"? We won't answer this question in detail here but provide some examples and intuition for their nature. The above mentioned consistency conditions can be loosely classified as a) reflecting the dynamics of the underlying process and b) the support of its distribution. Conditions of type a) are affine relations that the moments of process have to satisfy. To derive these conditions, note that the (extended) infinitesimal generator
$$
\mathcal{A} : w(t,z) \to \lim_{h\to 0^+} \frac{\mathbb{E_z[w(t + h,x(t + h))]} - w(t,z)}{h} = \frac{ \partial w(t,z) }{\partial t } + f(z)^\top \nabla_z w(t,z) + \text{Tr}\left(gg^\top(z) \nabla_z^2 w(t,z) \right) + \sum_{i=1}^{n_R} a_i(z) w(t, h_i(z))
$$
maps polynomials to polynomials under the assumption of a polynomial jump diffusion process. The moments of a polynomial jump-diffusion process accordingly follow linear, albeit generally underdetermined, dynamics:
$$
\frac{d}{dt}\mathbb{E}\left[\prod_{k=1}^n x_k^{i_k}(t) \right] = \mathbb{E}\left[ \mathcal{A}\prod_{k=1}^n x_k^{i_k}(t) \right] \iff \frac{dy_i}{dt}(t) = a_i^\top \mathbf{y}_q(t)
$$
From this point we could in principle derive more tractable conditions; however, since this requires introduction of more technical terms and we believe it does not contribute to building better intuition, we will refer the interested reader to the references listed below instead of going through the construction here.

Conditions of type b) impose positive semidefiniteness of certain moment matrices. To understand why such conditions are in fact somehow natural, consider a one-dimensional process $x(t)$ and a vector of the monomial basis $b(x) = [1, x, x^2, \dots, x^d]$. Then clearly,  the moment matrix
$$
\mathbb{E}\left[ b(x(t)) b(x(t))^\top \right] = \begin{bmatrix} 1 & y_1(t) & y_2(t)&\cdots & y_d(t) \\
															    y_1(t) & y_2(t) & y_3(t) & \cdots & y_{d+1}(t) \\
															    y_2(t) & y_3(t)  & \ddots & & y_{d+2}(t) \\
															    \vdots & \vdots & &  &\vdots \\
															    y_d(t) & y_{d+1}(t) & y_{d+2}(t) & \cdots & y_{2d}(t) 
															    \end{bmatrix}
$$
must be positive semidefinite as the left-hand-side is. This argument generalizes immediately to the multivariate case. The condition,
$$
\mathbb{E}\left[ b(x(t)) b(x(t))^\top \right]= \int_X b(x) b(x)^\top \, dP(x,t) \succeq 0,
$$
can be viewed as reflecting non-negativity of the probability measure $P(\cdot,t)$ describing the distribution of the process state at time $t$. With this intuition in mind, it follows further that for any polynomial $p$ which is non-negative on the state space $X$, the condition
$$
\mathbb{E}[p(x(t))b(x(t)) b(x(t))^\top] \succeq 0
$$
must hold, reflecting the support of the probability distribution $P(\cdot,t)$ on the state space $X$. Further observe that conditions of this form translate directly into [Linear Matrix Inequalities](https://en.wikipedia.org/wiki/Linear_matrix_inequality) on the moments of the process, suggesting why the resulting problems can be tackled via SDP. 

For more details and technicalities on moment bounding schemes, please consult the references below.

# Bounds on Stationary Moments of Markov Processes

```@docs
stationary_polynomial(MP::MarkovProcess, v::APL, d::Int, solver)
stationary_mean(MP::MarkovProcess, v::APL, d::Int, solver)
stationary_mean(rn::ReactionSystem, S0::Dict, S, d::Int, solver,
                scales = Dict(s => 1 for s in speceies(rn));
                auto_scaling = false)
stationary_variance(MP::MarkovProcess, v::APL, d::Int, solver)
stationary_variance(rn::ReactionSystem, S0, x, d::Int, solver,
                    scales = Dict(s => 1 for s in speceies(rn));
                    auto_scaling = false)
stationary_covariance_ellipsoid(MP::MarkovProcess, v::Vector{<:APL}, d::Int, solver)
stationary_covariance_ellipsoid(rn::ReactionSystem, S0::Dict, S::AbstractVector, d::Int, solver,
                                scales = Dict(s => 1 for s in speceies(rn));
                                auto_scaling = false)
```

# Bounds on Transient Moments of Markov Processes
```@docs
transient_polynomial(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)
transient_mean(MP::MarkovProcess, μ0::Dict, x::APL, d::Int, trange::AbstractVector{<:Real}, solver)
transient_mean(rn::ReactionSystem, S0::Dict, S, d::Int, trange::AbstractVector{<:Number}, solver,
            scales = Dict(s => 1 for s in speceies(rn));
            auto_scaling = false)
transient_variance(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)
transient_variance(rn::ReactionSystem, S0::Dict, S, d::Int, trange::AbstractVector{<:Real}, solver,
            scales = Dict(s => 1 for s in speceies(rn));
            auto_scaling = false)
transient_covariance_ellipsoid(MP::MarkovProcess, μ0::Dict, v::Vector{APL}, d::Int, trange::AbstractVector{<:Real}, solver)
transient_covariance_ellipsoid(rn::ReactionSystem, S0::Dict, S::AbstractVector, d::Int, trange::AbstractVector{<:Real}, solver,
            scales = Dict(s => 1 for s in speceies(rn));
            auto_scaling = false)
```

# Bounds on Stochastic Optimal Control Problems
```@docs
optimal_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver)
```
