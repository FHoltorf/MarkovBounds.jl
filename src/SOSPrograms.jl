## Stationary formulations

"""
	stationary_pop(MP::MarkovProcess, v::APL, d::Int, solver)

returns SOS program of degree d for compuitation of a *lower* bound
on the expecation of a polynomial observable v(x) at steady state of the
Markov process MP.
"""
function stationary_pop(MP::MarkovProcess, v::APL, d::Int, solver)
    model = SOSModel(solver)
    @variable(model, w, Poly(monomials(MP.x, 0:d)))
    @variable(model, s)
    @constraint(model, v + inf_generator(MP, w) >= s, domain = MP.X)
    @objective(model, Max, s)
    return model
end

"""
	stationary_polynomial(MP::MarkovProcess, v::APL, d::Int, solver)

returns a *lower* bound on the expecation of a polynomial observables v(x)
at steady state of the Markov process MP. The bound is computed based on an
SOS program over a polynomial of degree at most d; the bounds can be
tightened by increasing d. The program is solved with solver.
"""
function stationary_polynomial(MP::MarkovProcess, v::APL, d::Int, solver)
    model = stationary_pop(MP, v, d, solver)
    optimize!(model)
    return objective_value(model), termination_status(model), MOI.get(model, MOI.SolveTime())
end


"""
	stationary_mean(MP::MarkovProcess, v::APL, d::Int, solver)

returns *lower* and *upper* bound on the observable v(x) at steady state of
the Markov process MP. Bounds are computed based on SOS programs over a
polynomial of degree at most d; the bounds can be tightened by increasing
d. The program is solved with solver.
"""
function stationary_mean(MP::MarkovProcess, v::APL, d::Int, solver)
    lb, lb_stat, lb_time = stationary_polynomial(MP, v, d, solver)
    ub, ub_stat, ub_time = stationary_polynomial(MP, -v, d, solver)
    return [lb, - ub], [lb_stat,ub_stat], [lb_time, ub_time]
end

stationary_mean(RP::ReactionProcess, S, d::Int, solver) = stationary_mean(RP.JumpProcess, RP.species_to_state[S], d, solver)

"""
	stationary_mean(rn::ReactionSystem, S0::Dict, S, d::Int, solver,
			scales = Dict(s => 1 for s in species(rn));
			auto_scaling = false)

returns *lower* and *upper* bound on the mean of species S of the reaction
network rn with initial condition S0 (for all species!). The bound is based
on an SOS program of order d solved via solver; the bounds can be tightened
by increasing d.

For numerical stability, it is recommended to provide scales of the expected
magnitude of molecular counts for the different species at steady state.
If the system is *closed* it is also possible to enable auto_scaling which will
find the maximum molecular counts for each species under stoichiometry
constraints (via LP).

If the initial condition of the reaction network under investigation is
unknown or irrelevant, simply call

	stationary_mean(rn::ReactionSystem, S, d::Int, solver,
			scales = Dict(s => 1 for s in species(rn))).
"""
function stationary_mean(rn::ReactionSystem, S0::Dict, S, d::Int, solver,
	 					 scales = Dict(s => 1 for s in species(rn));
						 auto_scaling = false)
	RP, S0 = setup_reaction_process(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver)
 	return stationary_mean(RP.JumpProcess, RP.species_to_state[S], d, solver)
end
function stationary_mean(rn::ReactionSystem, S, d::Int, solver,
					     scales = Dict(s => 1 for s in species(rn)))
	RP = setup_reaction_process(rn, scales = scales)
	return stationary_mean(RP.JumpProcess, RP.species_to_state[S], d, solver)
end


"""
	stationary_variance(MP::MarkovProcess, v::APL, d::Int, solver)

returns SOS program of degree d for computation of an *upper* bound on the
variance of a polynomial observables v at steady state of the Markov process
MP.
"""
function stationary_variance(MP::MarkovProcess, v::APL, d::Int, solver)
    model = SOSModel(solver)
    @variable(model, w, Poly(monomials(MP.x, 0:d)))
    @variable(model, s)
    @variable(model, S[1:2])
    @constraint(model, [1 S[1]; S[1] S[2]] in PSDCone())
    @constraint(model, -(v^2 + 2*S[1]*v) + inf_generator(MP, w) >= s, domain = MP.X)
    @objective(model, Max, s-S[2])
    optimize!(model)
    return -objective_value(model), termination_status(model), MOI.get(model, MOI.SolveTime())
end

stationary_variance(RP::ReactionProcess, x, d::Int, solver) = stationary_variance(RP.JumpProcess, RP.species_to_state[x], d, solver)

"""
	stationary_variance(rn::ReactionSystem, S0, x, d::Int, solver,
			    scales = Dict(s => 1 for s in species(rn));
			    auto_scaling = false)

returns *upper* bound on the variance of species S of the reaction
network rn with initial condition S0 (for all species!). The bound is based
on an SOS program of degree d solved via solver; the bound can be tightened
by increasing d.

For numerical stability, it is recommended to provide scales of the expected
magnitude of molecular counts for the different species at steady state. If
the system is *closed* it is also possible to enable auto_scaling which will
find the maximum molecular counts for each species under stoichiometry
constraints (via LP).

If the initial condition of the reaction network under investigation is
unknown or irrelevant, simply call

	stationary_variance(rn::ReactionSystem, S, d::Int, solver,
			    scales = Dict(s => 1 for s in species(rn)))
"""
function stationary_variance(rn::ReactionSystem, S0, S, d::Int, solver,
	 						 scales = Dict(s => 1 for s in species(rn));
							 auto_scaling = false)
	RP, S0 = setup_reaction_process(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver)
 	return stationary_variance(RP.JumpProcess, RP.species_to_state[S], d, solver)
end
function stationary_variance(rn::ReactionSystem, S, d::Int, solver,
							 scales = Dict(s => 1 for s in species(rn)))
	RP = setup_reaction_process(rn, scales = scales)
	return stationary_variance(RP.JumpProcess, RP.species_to_state[S], d, solver)
end

"""
	stationary_covariance_ellipsoid(MP::MarkovProcess, v::Vector{<:APL}, d::Int, solver)

returns an *upper* on the volume of the covariance ellipsoid of a vector of
polynomial observables v(x), i.e., det(𝔼(v(x)v(x)ᵀ)), at steady state of the
Markov process MP.
The bounds are computed via an SOS program of degree d, hence can be tightened
by increasing d. This computation requires a solver that can deal with
exponential cone constraints.
"""
function stationary_covariance_ellipsoid(MP::MarkovProcess, v::Vector{<:APL}, d::Int, solver)
    n = length(v)
    model = SOSModel(solver)
    @variable(model, w, Poly(monomials(MP.x, 0:d)))
    @variable(model, s)
    @variable(model, S[1:n+1, 1:n+1], PSD)
    @variable(model, U[1:2*n, 1:2*n], PSD)
    @variable(model, r[1:n])
    @variable(model, q[1:n])
    @constraint(model, -(v'*S[1:n, 1:n]*v + 2*S[end,1:n]'*v) + inf_generator(MP, w) >= s, domain = MP.X)
    @constraint(model, S[1:n,1:n] .== U[1:n, 1:n])
    @constraint(model, [i in 1:n, j in 1:i], -2*U[n+i,j] - (i == j ? U[n+i,n+i] - r[i] : 0) == 0)
    @constraint(model, [i in 1:n], [-1, q[i], r[i]] in MOI.DualExponentialCone())
    @objective(model, Max, s-S[n+1,n+1]-sum(q))
    optimize!(model)
    return exp(-objective_value(model)), termination_status(model), MOI.get(model, MOI.SolveTime())
end

stationary_covariance_ellipsoid(RP::ReactionProcess, v, d::Int, solver) = stationary_variance(RP.JumpProcess, [RP.species_to_state[x] for x in v], d, solver)

"""
	stationary_covariance_ellipsoid(rn::ReactionSystem, S0::Dict, S::AbstractVector, d::Int, solver,
					scales = Dict(s => 1 for s in species(rn));
					auto_scaling = false)

returns an *upper* on the volume of the covariance ellipsoid of any subset S
of the chemical species in the reaction network rn, i.e., det(𝔼(SSᵀ)),
at steady state of the associated jump process. The reaction network is assumed
to have the deterministic initial state S0 (all species must be included here!).
The bounds are computed via an SOS program of degree d, hence can be tightened
by increasing d. This computation requires a solver that can deal with
exponential cone constraints.

For numerical stability, it is recommended to provide scales of the expected
magnitude of molecular counts for the different species at steady state. If
the system is *closed* it is also possible to enable auto_scaling which will
find the maximum molecular counts for each species under stoichiometry
constraints (via LP).

If the initial condition of the reaction network under investigation is
unknown or irrelevant, simply call

	stationary_covariance_ellipsoid(rn::ReactionSystem, S, d::Int, solver,
					scales = Dict(s => 1 for s in species(rn)))
"""
function stationary_covariance_ellipsoid(rn::ReactionSystem, S0::Dict, S::AbstractVector, d::Int, solver,
	 									 scales = Dict(s => 1 for s in species(rn));
										 auto_scaling = false)
	RP, x0 = setup_reaction_process(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver)
 	return stationary_covariance_ellipsoid(RP.JumpProcess, [RP.species_to_state[x] for x in S], d, solver)
end

function stationary_covariance_ellipsoid(rn::ReactionSystem, S, d::Int, solver,
	 									 scales = Dict(s => 1 for s in species(rn)))
	RP = setup_reaction_process(rn, scales = scales)
	return stationary_covariance_ellipsoid(RP.JumpProcess, [RP.species_to_state[x] for x in S], d, solver)
end

function stationary_probability_mass(MP::MarkovProcess, X::BasicSemialgebraicSet, d::Int, solver)
    ## TODO: Need good way to compute complement as union of several sets
end

## transient formulations

"""
	transient_pop(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)

returns SOS program of degree d for computing a *lower* bound on 𝔼[v(x(T))]
where v is a polynomial and x(T) the state of the Markov process MP at time T.
μ0 encodes the distribution of the initial state of the process in terms of its moments;
specifically, it maps monomials to the respective moments of the initial distribution.
trange is an *ordered* collection of time points used to discretize the time
horizon [0,T], i.e., trange[end] = T. Populating trange and increasing d has a
tightening effect on the bound furnished by the assembled SOS program.
"""
function transient_pop(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)
	if trange[1] == 0
		trange = trange[2:end]
	end
    nT = length(trange)
    @polyvar(t)
    T = @set(t >= 0 && t <= 1)
    XT = intersect(MP.X, T)
    Δt = [(i == 1 ? trange[i] : trange[i] - trange[i-1]) for i in 1:nT]
    model = SOSModel(solver)
    @variable(model, w[1:nT], Poly(monomials(sort([MP.x...,t], rev = true), 0:d)))
    @constraint(model, [i in 1:nT], extended_inf_generator(MP, w[i], t, scale=Δt[i]) >= 0, domain = XT)
    @constraint(model, [i in 1:nT-1], subs(w[i+1], t => 0) - subs(w[i], t => 1) >= 0, domain = MP.X)
    @constraint(model, v - subs(w[nT], t=>1) >= 0, domain = MP.X)
    @objective(model, Max, expectation(polynomial(subs(w[1], t=>0)), μ0))
    return model
end


"""
	transient_polynomial(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)

returns a *lower* bound on 𝔼[v(x(T))] where v is a polynomial and x(T) the state
of the Markov process MP at time T. μ0 encodes the
distribution of the initial state of the process in terms of its moments; specifically,
it maps monomials to the respective moments of the initial distribution.
trange is an *ordered* collection of time points used to discretize the time
horizon [0,T], i.e., trange[end] = T. Populating trange and increasing d improves
the computed bound.
"""
function transient_polynomial(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)
    model = transient_pop(MP, μ0, v, d, trange, solver)
    optimize!(model)
    return objective_value(model), termination_status(model), MOI.get(model, MOI.SolveTime())
end

"""
	transient_mean(MP::MarkovProcess, μ0::Dict, x::APL, d::Int, trange::AbstractVector{<:Real}, solver)

returns a *lower* and *upper* bound on 𝔼[v(x(T))] where v is a polynomial and x(T) the state
of the Markov process MP at time T. μ0 encodes the
distribution of the initial state of the process in terms of its moments; specifically,
it maps monomials to the respective moments of the initial distribution.
trange is an *ordered* collection of time points used to discretize the time
horizon [0,T], i.e., trange[end] = T. Populating trange and increasing d improves
the computed bounds.
"""
function transient_mean(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)
    lb, lb_stat, lb_time = transient_polynomial(MP, μ0, v, d, trange, solver)
    ub, ub_stat, ub_time = transient_polynomial(MP, μ0, -v, d, trange, solver)
    return [lb, -ub], [lb_stat,ub_stat], [lb_time, ub_time]
end

"""
	transient_mean(rn::ReactionSystem, S0::Dict, S, d::Int, trange::AbstractVector{<:Number}, solver,
			scales = Dict(s => 1 for s in species(rn));
			auto_scaling = false)

returns a *lower* and *upper* bound on the mean of the molecular count of species
S in reaction network rn at time T. S0 refers to the *deterministic* initial state
of the reaction system (including all species!).
trange is an *ordered* collection of time points used to discretize the time
horizon [0,T], i.e., trange[end] = T. Populating trange and increasing d improves
the computed bounds.

For numerical stability, it is recommended to provide scales of the expected
magnitude of molecular counts for the different species at steady state. If
the system is *closed* it is also possible to enable auto_scaling which will
find the maximum molecular counts for each species under stoichiometry
constraints (via LP).
"""
function transient_mean(rn::ReactionSystem, S0::Dict, S, d::Int, trange::AbstractVector{<:Number}, solver,
		  				scales = Dict(s => 1 for s in species(rn));
						auto_scaling = false)
	RP, S0 = setup_reaction_process(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver)
	μ0 = init_moments(RP.JumpProcess.x, S0, d + maximum(maxdegree.(RP.JumpProcess.a)))
 	return transient_mean(RP.JumpProcess, μ0, RP.species_to_state[S], d, trange, solver)
end


"""
	transient_variance(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)

returns an *upper* bound on 𝔼[v(x(T))²] - 𝔼[v(x(T))]² where v is a polynomial
and x(T) the state of the Markov process MP at time T.
"""
function transient_variance(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)
	if trange[1] == 0
		trange = trange[2:end]
	end
	nT = length(trange)
	@polyvar(t)
    Δt = [(i == 1 ? trange[i] : trange[i] - trange[i-1]) for i in 1:nT]
	T = @set(t >= 0 && t <= 1)
    XT = intersect(MP.X, T)

	model = SOSModel(solver)
	@variable(model, w[1:nT], Poly(monomials(sort([MP.x...,t], rev = true), 0:d)))
	@variable(model, S[1:2])
	@constraint(model, [i in 1:nT], extended_inf_generator(MP, w[i], t, scale = Δt[i]) >= 0, domain = XT)
	@constraint(model, [i in 1:nT-1], subs(w[i+1], t => 0) - subs(w[i], t => 1) >= 0, domain = MP.X)
	@constraint(model, - (v^2 + 2*S[1]*v) - subs(w[nT], t => 1) >= 0, domain = MP.X)
	@constraint(model, [1 S[1]; S[1] S[2]] in PSDCone())
	@objective(model, Max, expectation(polynomial(subs(w[1], t => 0)), μ0) - S[2])
	optimize!(model)
	return -objective_value(model), termination_status(model), MOI.get(model, MOI.SolveTime())
end

"""
	transient_variance(rn::ReactionSystem, S0::Dict, S, d::Int, trange::AbstractVector{<:Real}, solver,
			    scales = Dict(s => 1 for s in species(rn));
				auto_scaling = false)

returns an *upper* bound on the variance of species S in the reaction network rn
at time T. S0 refers to the *deterministic* initial state
of the reaction system (including all species!).
trange is an *ordered* collection of time points used to discretize the time
horizon [0,T], i.e., trange[end] = T. Populating trange and increasing d improves
the computed bounds.

For numerical stability, it is recommended to provide scales of the expected
magnitude of molecular counts for the different species at steady state. If
the system is *closed* it is also possible to enable auto_scaling which will
find the maximum molecular counts for each species under stoichiometry
constraints (via LP).
"""
function transient_variance(rn::ReactionSystem, S0::Dict, S, d::Int, trange::AbstractVector{<:Real}, solver,
						    scales = Dict(s => 1 for s in species(rn));
							auto_scaling = false)
	RP, S0 = setup_reaction_process(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver)
	μ0 = init_moments(RP.JumpProcess.x, S0, d + maximum(maxdegree.(RP.JumpProcess.a)))
 	return transient_variance(RP.JumpProcess, μ0, RP.species_to_state[S], d, trange, solver)
end

"""
	transient_covariance_ellipsoid(MP::MarkovProcess, μ0::Dict, v::Vector{APL}, d::Int, trange::AbstractVector{<:Real}, solver)

returns an *upper* bound on the volume of the covariance ellipsoid det(𝔼(v(x(T))v(x(T))ᵀ)),
where v is a polynomial and x(T) the state of the Markov process MP at time T.
"""
function transient_covariance_ellipsoid(MP::MarkovProcess, μ0::Dict, v::Vector{<:APL}, d::Int, trange::AbstractVector{<:Real}, solver)
	if trange[1] == 0
		trange = trange[2:end]
	end
	nT = length(trange)
	n = length(v)
	@polyvar(t)
    Δt = [(i == 1 ? trange[i] : trange[i] - trange[i-1]) for i in 1:nT]
	T = @set(t >= 0 && t <= 1)
    XT = intersect(MP.X, T)

	model = SOSModel(solver)
	@variable(model, w[1:nT], Poly(monomials(sort([MP.x...,t], rev = true), 0:d)))
	@variable(model, S[1:n+1,1:n+1], PSD)
	@variable(model, U[1:2n,1:2n], PSD)
	@variable(model, r[1:n])
    @variable(model, q[1:n])
	@constraint(model, [i in 1:nT], extended_inf_generator(MP, w[i], t, scale = Δt[i]) >= 0, domain = XT)
	@constraint(model, [i in 1:nT-1], subs(w[i+1], t => 0) - subs(w[i], t => 1) >= 0, domain = MP.X)
	@constraint(model, -(v'*S[1:n, 1:n]*v + 2*S[end,1:n]'*v) - subs(w[nT], t=>1) >= 0, domain = MP.X)
    @constraint(model, S[1:n,1:n] .== U[1:n, 1:n])
    @constraint(model, [i in 1:n, j in 1:i], -2*U[n+i,j] - (i == j ? U[n+i,n+i] - r[i] : 0) == 0)
    @constraint(model, [i in 1:n], [-1, q[i], r[i]] in MOI.DualExponentialCone())
	@objective(model, Max, expectation(polynomial(subs(w[1], t => 0)), μ0) - S[n+1,n+1] - sum(q))
	optimize!(model)
	return exp(-objective_value(model)), termination_status(model), MOI.get(model, MOI.SolveTime())
end

"""
	transient_covariance_ellipsoid(rn::ReactionSystem, S0::Dict, S::AbstractVector, d::Int, trange::AbstractVector{<:Real}, solver,
							scales = Dict(s => 1 for s in species(rn));
							auto_scaling = false)

returns an *upper* bound on the volume of the covariance ellipsoid associated with any
collection of chemical species in the reaction network rn at time T.
S0 refers to the *deterministic* initial state of the reaction system
(including all species!). trange is an *ordered* collection of time points used to discretize the time
horizon [0,T], i.e., trange[end] = T. Populating trange and increasing d improves
the computed bounds.

For numerical stability, it is recommended to provide scales of the expected
magnitude of molecular counts for the different species at steady state. If
the system is *closed* it is also possible to enable auto_scaling which will
find the maximum molecular counts for each species under stoichiometry
constraints (via LP).
"""
function transient_covariance_ellipsoid(rn::ReactionSystem, S0::Dict, S::AbstractVector, d::Int, trange::AbstractVector{<:Real}, solver,
										scales = Dict(s => 1 for s in species(rn));
										auto_scaling = false)
	RP, x0 = setup_reaction_process(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver)
	μ0 = init_moments(RP.JumpProcess.x, x0, d + maximum(maxdegree.(RP.JumpProcess.a)))
 	return transient_covariance_ellipsoid(RP.JumpProcess, μ0, [RP.species_to_state[x] for x in S], d, trange, solver)
end


## Optimal Control
"""
	optimal_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver)

returns a *lower* bound on the objective value of the (stochastic) optimal
control problem specified by CP. μ0 encodes information about the distribution of
the initial state of the process; specifically, μ0 maps a given monomial to the
corresponding moment of the initial distribution. trange refers to an *ordered*
set of time points discretizing the control horizon. trange[end] should coincide
with the end of the control horizon, i.e., trange[end] = Inf in case of an infinite
horizon problem. The bound is computed via a SOS program of degree d solved with
an appropriate method given by solver.

The bound can be tightened by populating trange or increasing d.
"""
function optimal_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver)
	if isinf(CP.T)
		model = infinite_horizon_control(CP, μ0, d, trange, solver)
	else
		model = finite_horizon_control(CP, μ0, d, trange, solver)
	end
	optimize!(model)
	return  objective_value(model), termination_status(model), MOI.get(model, MOI.SolveTime()), model
end

## Finite horizon control problems
function finite_horizon_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver)
	if isa(CP.Objective, LagrangeMayer)
		model = finite_horizon_LM(CP, μ0, d, trange, solver)
	elseif isa(CP.Objective, ExitProbability)
		model = finite_horizon_EP(CP, μ0, d, trange, solver)
	elseif isa(CP.Objective, TerminalSetProbability)
		model = finite_horizon_TP(CP, μ0, d, trange, solver)
	else
		@error "Objective function type not supported"
	end
	return model
end

function finite_horizon_LM(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver)
	if trange[1] == 0
		trange = trange[2:end]
	end
	nT = length(trange)
	MP = CP.MP
	Δt = [(i == 1 ? trange[i] : trange[i] - trange[i-1]) for i in 1:nT]
	T = @set(CP.t >= 0 && CP.t <= 1)
	XTU = intersect(MP.X,T,CP.U)

	model = SOSModel(solver)
	@variable(model, w[1:nT], Poly(monomials(sort([MP.x...,CP.t],rev=true), 0:d)))
	@constraint(model, [i in 1:nT], CP.Objective.l + extended_inf_generator(MP, w[i], CP.t, scale = Δt[i]) >= 0, domain = XTU)
	@constraint(model, [i in 1:nT-1], subs(w[i+1], CP.t => 0) - subs(w[i], CP.t => 1) >= 0, domain = MP.X)
	@constraint(model, CP.Objective.m - subs(w[nT], CP.t=>1) >= 0, domain = MP.X)
	obj = expectation(polynomial(subs(w[1], CP.t => 0)), μ0)
	if !isempty(CP.PathChanceConstraints)
		for PC in CP.PathChanceConstraints
			∂X = ∂(PC.X)
			s_p = @variable(model)
			v = @variable(model, [1:nT], Poly(monomials(sort([MP.x..., CP.t],rev=true), 0:d)))
			@constraint(model, [i in 1:nT], extended_inf_generator(MP, v[i], CP.t, scale = Δt[i]) >= 0, domain = intersect(PC.X, T, CP.U))
			@constraint(model, [i in 1:nT-1], subs(v[i+1], CP.t => 0) - subs(v[i], CP.t => 1) >= 0, domain = intersect(PC.X, T))
			@constraint(model, [i in 1:nT, k in 1:length(∂X)], v[i] >= 0, domain = ∂X[k])
			@constraint(model, -s_p - subs(v[nT], CP.t => 1) >= 0, domain = PC.X)
			@constraint(model, s_p >= 0)
			obj += s_p*(1-PC.α)
		end
	end
	if !isempty(CP.TerminalChanceConstraints)
		for TC in CP.TerminalChanceConstraints
			s_p = @variable(model)
			@constraint(model, s_p >= 0)
			@constraint(model, - s_p + CP.Objective.m - subs(w[nT], CP.t => 1) >= 0, domain = TC.X)
			obj += s_p*(1-TC.α)
		end
	end
	@objective(model, Max, obj)
	return model
end

function finite_horizon_EP(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver)
	if trange[1] == 0
		trange = trange[2:end]
	end
	nT = length(trange)
	MP = CP.MP
	Δt = [(i == 1 ? trange[i] : trange[i] - trange[i-1]) for i in 1:nT]
	T = @set(CP.t >= 0 && CP.t <= 1)
	XTU = intersect(CP.Objective.X,MP.X,T,CP.U)

	model = SOSModel(solver)
	∂X = ∂(CP.Objective.X)
	@variable(model, s >= 0)
	@variable(model, w[1:nT], Poly(monomials(sort([MP.x...,CP.t],rev=true), 0:d)))
	@constraint(model, [i in 1:nT], extended_inf_generator(MP, w[i], CP.t, scale = Δt[i]) >= 0, domain = XTU)
	@constraint(model, [i in 1:nT-1], subs(w[i+1], CP.t => 0) - subs(w[i], CP.t => 1) >= 0, domain = intersect(CP.Objective.X, MP.X))
	@constraint(model, [i in 1:nT, k in 1:length(∂X)], w[i] >= 0, domain = intersect(∂X[k], MP.X, T))
	@constraint(model, - 1 - subs(w[nT], CP.t => 1) >= 0, domain = CP.Objective.X)
	obj = expectation(polynomial(subs(w[1], CP.t => 0)), μ0)
	if !isempty(CP.PathChanceConstraints)
		for PC in CP.PathChanceConstraints
			∂X = ∂(PC.X)
			s_p = @variable(model)
			v = @variable(model, [1:nT], Poly(monomials(sort([MP.x..., CP.t],rev=true), 0:d)))
			@constraint(model, [i in 1:nT], extended_inf_generator(MP, v[i], CP.t, scale = Δt[i]) >= 0, domain = intersect(PC.X, T, CP.U))
			@constraint(model, [i in 1:nT-1], subs(v[i+1], CP.t => 0) - subs(v[i], CP.t => 1) >= 0, domain = PC.X)
			@constraint(model, [i in 1:nT, k in 1:length(∂X)], v[i] >= 0, domain = intersect(∂X[k], T))
			@constraint(model, -s_p - subs(v[nT], CP.t => 1) >= 0, domain = PC.X)
			@constraint(model, s_p >= 0)
			obj += s_p*(1-PC.α)
		end
	end
	if !isempty(CP.TerminalChanceConstraints)
		@error "Terminal chance constraints currently not supported for Exit Probability problems!"
	end
	@objective(model, Max, obj)
	return model
end

function finite_horizon_TP(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver)
	if trange[1] == 0
		trange = trange[2:end]
	end
	nT = length(trange)
	MP = CP.MP
	Δt = [(i == 1 ? trange[i] : trange[i] - trange[i-1]) for i in 1:nT]
	T = @set(CP.t >= 0 && CP.t <= 1)
	XTU = intersect(MP.X,T,CP.U)

	model = SOSModel(solver)
	@variable(model, w[1:nT], Poly(monomials(sort([MP.x...,CP.t],rev=true), 0:d)))
	@constraint(model, [i in 1:nT], extended_inf_generator(MP, w[i], CP.t, scale = Δt[i]) >= 0, domain = XTU)
	@constraint(model, [i in 1:nT-1], subs(w[i+1], CP.t => 0) - subs(w[i], CP.t => 1) >= 0, domain = MP.X)
	@constraint(model, - subs(w[nT], CP.t=>1) >= 0, domain = MP.X)
	@constraint(model, - 1 - subs(w[nT], CP.t=>1) >= 0, domain = CP.Objective.X)
	obj = expectation(polynomial(subs(w[1], CP.t => 0)), μ0)
	if !isempty(CP.PathChanceConstraints)
		for PC in CP.PathChanceConstraints
			∂X = ∂(PC.X)
			s_p = @variable(model)
			v = @variable(model, [1:nT], Poly(monomials(sort([MP.x..., CP.t],rev=true), 0:d)))
			@constraint(model, [i in 1:nT], extended_inf_generator(MP, v[i], CP.t, scale = Δt[i]) >= 0, domain = intersect(PC.X, T, CP.U))
			@constraint(model, [i in 1:nT-1], subs(v[i+1], CP.t => 0) - subs(v[i], CP.t => 1) >= 0, domain = intersect(PC.X, T))
			@constraint(model, [i in 1:nT, k in 1:length(∂X)], v[i] >= 0, domain = ∂X[k])
			@constraint(model, -s_p - subs(v[nT], CP.t => 1) >= 0, domain = PC.X)
			@constraint(model, s_p >= 0)
			obj += s_p*(1-PC.α)
		end
	end
	if !isempty(CP.TerminalChanceConstraints)
		for TC in CP.TerminalChanceConstraints
			s_p = @variable(model)
			@constraint(model, s_p >= 0)
			@constraint(model, - s_p + CP.Objective.m - subs(w[nT], CP.t => 1) >= 0, domain = TC.X)
			obj += s_p*(1-TC.α)
		end
	end
	@objective(model, Max, obj)
	return model
end

## Discounted inf horizon control problems
function infinite_horizon_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver)
	@assert isa(CP.Objective, LagrangeMayer) "Objective function type not supported"
	@assert CP.Objective.m == 0 "Mayer term must be 0 for infinite horizon problem"
	return infinite_horizon_LM(CP, μ0, d, trange, solver)
end

function infinite_horizon_LM(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver)
	if trange[1] == 0
		trange = trange[2:end]
	end
	nT = length(trange)
	MP = CP.MP
	ρ = CP.discount_factor
	Δt = [(i == 1 ? trange[i] : trange[i] - trange[i-1]) for i in 1:nT]
	T = @set(CP.t >= 0 && CP.t <= 1)
	XTU = intersect(MP.X,T,CP.U)

	model = SOSModel(solver)
	@variable(model, w[1:nT-1], Poly(monomials(sort([MP.x...,CP.t],rev=true), 0:d)))
	@variable(model, w[nT], Poly(monomials(MP.x, 0:d)))
	@constraint(model, [i in 1:nT-1], CP.Objective.l + extended_inf_generator(MP, w[i], CP.t, scale = Δt[i]) - ρ*w[i] >= 0, domain = XTU)
	@constraint(model, CP.Objective.l + inf_generator(MP, w[nT]) >= 0, domain = intersect(MP.X, CP.U))
	@constraint(model, - subs(w[nT]) >= 0, domain = MP.X)
	obj = expectation(polynomial(subs(w[1], CP.t => 0)), μ0)
	if !isempty(CP.PathChanceConstraints)
		@warn "Path chance constraints are omitted"
	end
	if !isempty(CP.TerminalChanceConstraints)
		@warn "Terminal chance constraints are omitted"
	end
	@objective(model, Max, obj)
	return model
end
