```@meta
CurrentModule = MarkovBounds
```

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
