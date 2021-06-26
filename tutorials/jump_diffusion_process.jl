using MarkovBoundsSOS, Plots, MosekTools
"""
    Beyond jump processes modeling stochastic chemical systems, a much wider range
    of stochastic processes admits computable bounds on stationary and transient
    moments and related statistics via convex optimization.
    Generally all jump-diffusion processes of the form
        dx = f(x) dt + g(x) dW + ∑ hᵢ(x) dn(aᵢ)
    can be analyzed as long as the data:
        f - drift coefficient
        ggᵀ - diffusion matrix
        hᵢ - jumps
        aᵢ - arrival rates
    are polynomials and the state space of the system is basic closed semialgebraic
    (i.e., described by finitely many polynomial inequalities) or can at least be
    reasonably well approximated by a basic closed semialgebraic set
    In the following, we show with an example how to set the bounding problems up
    and compute relevant statistics.

    Let us consider the simple CIR model for the dynamics of interest rates
        dx =  κ(θ-x) dt + σ √x dW
    As can be seen from the model this is a pure diffusion process. Such a
    process is defined as follows:
"""
κ, θ, σ = 0.15, 0.03, 0.05 # model parameters
@polyvar(x) # states (interest rate)
f = κ*(θ-x) # drift coefficient
g = σ^2*x   # diffusion matrix (outer product )
X = @set(x >= 0) # support/state space

DP = DiffusionProcess(x, f, g, X)
"""
    With the diffusion process defined, we can bound interesting quantities such
    as the long term average or variance of the interest rate x simply by calling
"""
mean = stationary_mean(DP, x, 2, Mosek.Optimizer)
var = stationary_variance(DP, x, 2, Mosek.Optimizer)
"""
    similarly, bounds along a trajectory can be evaluated with ease.
    In order to showcase how also jump-diffusion are dealt with, let us assume
    that the interest rate drops to half of its value at random times characterized
    by a Poisson process with rate a(x) = 0.035 x. To define this process, we can
    simply define the jump component separately as a jump process:
"""
a = 0.035*x # arrival rate
h = x/2 # jump (interest rate jumps to half its value)
JP = JumpProcess(x, a, h, X)
"""
    The overall jump-diffusion process is then defined in terms of the jump and
    diffusion process:
"""
JDP = JumpDiffusionProcess(JP, DP)
"""
    Alternatively the process could also be defined in terms of the whole set of
    problem data:
        ``JDP = JumpDiffusionProcess(x, a, h, f, g, X)``

    Now we can for example study the evolution of means and variances of this process
    over time:
"""
Ts = [0.1, 0.25, 0.5, 1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10] # time points to probe mean and variance at
x0 = 0.01 # initial condition (deterministic - assumed)
order = 4 # relaxation order used
nT = 10 # number of time intervals used to discretize the time domain
μ0 = Dict(x^i => x0^i for i in 0:order+1) # moments of the initial distribution
var_bounds, mean_bounds = [], []
for T in Ts
    trange = range(0, T, length=nT + 1)
    mean = transient_mean(JDP, μ0, x, order, trange, Mosek.Optimizer)
    push!(mean_bounds, mean)
    var = transient_variance(JDP, μ0, x, order, trange, Mosek.Optimizer)
    push!(var_bounds, var)
end
p = plot(xlabel="time", ylabel="interest rate", legend=:bottomright)
plot!(p, Ts, [b[1].value for b in mean_bounds], color=:blue, label="mean lower bound")
plot!(p, Ts, [b[2].value for b in mean_bounds], color=:red, label="mean upper bound")
plot!(p, Ts, sqrt.([var.value for var in var_bounds]), color=:red, linestyle=:dash, label="std upper bound")
display(p)
