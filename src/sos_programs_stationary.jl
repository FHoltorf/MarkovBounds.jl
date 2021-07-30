export stationary_polynomial, stationary_mean, stationary_variance, stationary_covariance_ellipsoid

"""
	stationary_pop(MP::MarkovProcess, v::APL, d::Int, solver)

returns SOS program of degree `d` for compuitation of a **lower** bound
on the expecation of a polynomial observable ``v(x)`` at steady state of the
Markov process `MP`.
"""
function stationary_pop(MP::MarkovProcess, p::APL, order::Int, solver, P::Partition = trivial_partition(MP.X))
    model = SOSModel(solver)
    @variable(model, s)
    w = Dict(v => (props(P.graph, v)[:cell] isa Singleton ?
                   @variable(model) :
                   @variable(model, [1], Poly(monomials(MP.x, 0:order)))[1]) for v in vertices(P.graph))

    for v in vertices(P.graph)
        add_stationarity_constraints!(model, MP, v, P, props(P.graph, v)[:cell], w, s-p)
    end

    for e in edges(P.graph)
        add_coupling_constraints!(model, MP, e, P, w)
    end
    @objective(model, Max, s)
    return model, w
end

"""
	stationary_polynomial(MP::MarkovProcess, v::APL, d::Int, solver)

returns a **lower** bound on the expecation of a polynomial observables ``v(x)``
at steady state of the Markov process `MP`. The bound is computed based on an
SOS program over a polynomial of degree at most `d`; the bounds can be
tightened by increasing `d`. The program is solved with `solver`.
"""
function stationary_polynomial(MP::MarkovProcess, v::APL, d::Int, solver, P::Partition = trivial_partition(MP.X))
    model, w = stationary_pop(MP, v, d, solver, P)
    optimize!(model)
    return Bound(objective_value(model), model, P, Dict(key => value(w[key]) for key in keys(w)))
end

stationary_polynomial(MP::MarkovProcess, v::Num, d::Int, solver, P::Partition = trivial_partition(MP.X)) = stationary_polynomial(MP, polynomialize_expr(v, MP.poly_vars), d, solver, P)

"""
	stationary_mean(MP::MarkovProcess, v::APL, d::Int, solver)

returns **lower** and **upper** bound on the observable ``v(x)`` at steady state of
the Markov process `MP`. Bounds are computed based on SOS programs over a
polynomial of degree at most `d`; the bounds can be tightened by increasing
`d`. The program is solved with `solver`.
"""
function stationary_mean(MP::MarkovProcess, v::APL, d::Int, solver, P::Partition = trivial_partition(MP.X))
    lb = stationary_polynomial(MP, v, d, solver, P)
    ub = stationary_polynomial(MP, -v, d, solver, P)
	ub.value *= -1
    return lb, ub
end

stationary_mean(RP::ReactionProcess, S, d::Int, solver, P::Partition = trivial_partition(RP.JumpProcess.X)) = stationary_mean(RP.JumpProcess, RP.species_to_state[S], d, solver, P)
stationary_mean(LP::LangevinProcess, S, d::Int, solver, P::Partition = trivial_partition(LP.DiffusionProcess.X)) = stationary_mean(LP.DiffusionProcess, LP.species_to_state[S], d, solver, P)
stationary_mean(MP::MarkovProcess, v::Num, d::Int, solver, P::Partition = trivial_partition(MP.X)) = stationary_mean(MP, polynomialize_expr(v, MP.poly_vars), d, solver, P)

"""
	stationary_mean(rn::ReactionSystem, S0::Dict, S, d::Int, solver,
			scales = Dict(s => 1 for s in species(rn));
			auto_scaling = false)

returns **lower** and **upper** bound on the mean of species `S` of the reaction
network `rn` with initial condition `S0` (for all species!). The bound is based
on an SOS program of order `d` solved via `solver`; the bounds can be tightened
by increasing `d`.

For numerical stability, it is recommended to provide scales of the expected
magnitude of molecular counts for the different species at steady state.
If the system is **closed** it is also possible to enable auto_scaling which will
find the maximum molecular counts for each species under stoichiometry
constraints (via LP).

If the initial condition of the reaction network under investigation is
unknown or irrelevant, simply call

	stationary_mean(rn::ReactionSystem, S, d::Int, solver,
			scales = Dict(s => 1 for s in species(rn))).
"""
function stationary_mean(rn::ReactionSystem, S0::Dict, S, d::Int, solver,
	 					 scales = Dict(s => 1 for s in species(rn));
						 params::Dict = Dict(), auto_scaling = false)
	RP, S0 = reaction_process_setup(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver, params = params)
 	return stationary_mean(RP.JumpProcess, RP.species_to_state[S], d, solver)
end

function stationary_mean(rn::ReactionSystem, S, d::Int, solver,
					     scales = Dict(s => 1 for s in species(rn));
						 params::Dict = Dict())
	RP = reaction_process_setup(rn, scales = scales, params = params)
	return stationary_mean(RP.JumpProcess, RP.species_to_state[S], d, solver)
end


"""
	stationary_variance(MP::MarkovProcess, v::APL, d::Int, solver)

returns SOS program of degree `d` for computation of an **upper** bound on the
variance of a polynomial observables `v` at steady state of the Markov process
`MP`.
"""
function stationary_variance(MP::MarkovProcess, p::APL, d::Int, solver, P::Partition = trivial_partition(MP.X))
    model = SOSModel(solver)
	w = Dict(v => (props(P.graph, v)[:cell] isa Singleton ?
                   @variable(model) :
                   @variable(model, [1], Poly(monomials(MP.x, 0:d)))[1]) for v in vertices(P.graph))
    @variable(model, s)
    @variable(model, S[1:2])
    @constraint(model, [1 S[1]; S[1] S[2]] in PSDCone())

	for v in vertices(P.graph)
        add_stationarity_constraints!(model, MP, v, P, props(P.graph, v)[:cell], w, s + p^2 + 2*S[1]*p)
    end

	for e in edges(P.graph)
        add_coupling_constraints!(model, MP, e, P, w)
    end

    @objective(model, Max, s-S[2])
    optimize!(model)
    return Bound(-objective_value(model), model, P, Dict(key => value(w[key]) for key in keys(w)))
end

stationary_variance(RP::ReactionProcess, x, d::Int, solver, P::Partition = trivial_partition(RP.JumpProcess.X)) = stationary_variance(RP.JumpProcess, RP.species_to_state[x], d, solver, P)
stationary_variance(MP::MarkovProcess, x::Num, d::Int, solver, P::Partition = trivial_partition(MP.X)) = stationary_variance(MP, polynomialize_expr(x, MP.poly_vars), d, solver, P)

"""
	stationary_variance(rn::ReactionSystem, S0, x, d::Int, solver,
			    scales = Dict(s => 1 for s in species(rn));
			    auto_scaling = false)

returns **upper** bound on the variance of species `S` of the reaction
network rn with initial condition `S0` (for all species!). The bound is based
on an SOS program of degree `d` solved via `solver`; the bound can be tightened
by increasing `d`.

For numerical stability, it is recommended to provide scales of the expected
magnitude of molecular counts for the different species at steady state. If
the system is **closed** it is also possible to enable `auto_scaling` which will
find the maximum molecular counts for each species under stoichiometry
constraints (via LP).

If the initial condition of the reaction network under investigation is
unknown or irrelevant, simply call

	stationary_variance(rn::ReactionSystem, S, d::Int, solver,
			    scales = Dict(s => 1 for s in species(rn)))
"""
function stationary_variance(rn::ReactionSystem, S0, S, d::Int, solver,
	 						 scales = Dict(s => 1 for s in species(rn));
							 auto_scaling = false, params::Dict = Dict())
	RP, S0 = reaction_process_setup(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver, params = params)
 	return stationary_variance(RP.JumpProcess, RP.species_to_state[S], d, solver)
end

function stationary_variance(rn::ReactionSystem, S, d::Int, solver,
							 scales = Dict(s => 1 for s in species(rn));
							 params::Dict = Dict())
	RP = reaction_process_setup(rn, scales = scales, params = params)
	return stationary_variance(RP.JumpProcess, RP.species_to_state[S], d, solver)
end

@doc raw"""
	stationary_covariance_ellipsoid(MP::MarkovProcess, v::Vector{<:APL}, d::Int, solver)

returns an **upper** on the volume of the covariance ellipsoid of a vector of
polynomial observables ``v(x)``, i.e., ``\text{det}(\mathbb{E} [v(x)v(x)^\top] - \mathbb{E}[v(x)] \mathbb{E}[v(x)]^\top)``, at steady state of the
Markov process `MP`.
The bounds are computed via an SOS program of degree `d`, hence can be tightened
by increasing `d`. This computation requires a `solver` that can handle
exponential cone constraints.
"""
function stationary_covariance_ellipsoid(MP::MarkovProcess, v::Vector{<:APL}, d::Int, solver, P::Partition = trivial_partition(MP.X))
    n = length(v)
    model = SOSModel(solver)
	w = Dict(v => (props(P.graph, v)[:cell] isa Singleton ?
                   @variable(model) :
                   @variable(model, [1], Poly(monomials(MP.x, 0:d)))[1]) for v in vertices(P.graph))
    @variable(model, s)
    @variable(model, S[1:n+1, 1:n+1], PSD)
    @variable(model, U[1:2*n, 1:2*n], PSD)
    @variable(model, r[1:n])
    @variable(model, q[1:n])

	@constraint(model, S[1:n,1:n] .== U[1:n, 1:n])
    @constraint(model, [i in 1:n, j in 1:i], -2*U[n+i,j] - (i == j ? U[n+i,n+i] - r[i] : 0) == 0)
    @constraint(model, [i in 1:n], [-1, q[i], r[i]] in MOI.DualExponentialCone())

	for v in vertices(P.graph)
        flag = (v_target == v)
        add_stationarity_constraints!(model, MP, v, P, props(P.graph, v)[:cell], w, s + (v'*S[1:n, 1:n]*v + 2*S[end,1:n]'*v))
    end

	for e in edges(P.graph)
        add_coupling_constraints!(model, MP, e, P, w)
    end

    @objective(model, Max, s-S[n+1,n+1]-sum(q))
    optimize!(model)
    return Bound(exp(-objective_value(model)), model, P, Dict(key => value(w[key]) for key in keys(w)))
end

stationary_covariance_ellipsoid(RP::ReactionProcess, v::Vector, d::Int, solver, P::Partition = trivial_partition(RP.JumpProcess.X)) =
								stationary_variance(RP.JumpProcess, [RP.species_to_state[x] for x in v], d, solver, P)
stationary_covariance_ellipsoid(MP::MarkovProcess, v::Vector{Num}, d::Int, solver, P::Partition = trivial_partition(MP.X)) =
								stationary_variance(MP, polynomialize_expr(v, MP.poly_vars), d, solver, P)

@doc raw"""
	stationary_covariance_ellipsoid(rn::ReactionSystem, S0::Dict, S::AbstractVector, d::Int, solver,
					scales = Dict(s => 1 for s in species(rn));
					auto_scaling = false)

returns an **upper** on the volume of the covariance ellipsoid of any subset `S`
of the chemical species in the reaction network `rn`, i.e., ``\text{det}(\mathbb{E}[SS^\top] - \mathbb{E}[S] \mathbb{E}[S]^\top)``,
at steady state of the associated jump process. The reaction network is assumed
to have the deterministic initial state `S0` (all species must be included here!).
The bounds are computed via an SOS program of degree `d`, hence can be tightened
by increasing `d`. This computation requires a solver that can deal with
exponential cone constraints.

For numerical stability, it is recommended to provide scales of the expected
magnitude of molecular counts for the different species at steady state. If
the system is **closed** it is also possible to enable auto_scaling which will
find the maximum molecular counts for each species under stoichiometry
constraints (via LP).

If the initial condition of the reaction network under investigation is
unknown or irrelevant, simply call

	stationary_covariance_ellipsoid(rn::ReactionSystem, S, d::Int, solver,
					scales = Dict(s => 1 for s in species(rn)))
"""
function stationary_covariance_ellipsoid(rn::ReactionSystem, S0::Dict, S::AbstractVector, d::Int, solver,
	 									 scales = Dict(s => 1 for s in species(rn));
										 auto_scaling = false, params::Dict = Dict())
	RP, x0 = reaction_process_setup(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver, params = params)
 	return stationary_covariance_ellipsoid(RP.JumpProcess, [RP.species_to_state[x] for x in S], d, solver)
end

function stationary_covariance_ellipsoid(rn::ReactionSystem, S::Vector, d::Int, solver,
	 									 scales = Dict(s => 1 for s in species(rn));
										 params::Dict = Dict())
	RP = reaction_process_setup(rn, scales = scales, params = params)
	return stationary_covariance_ellipsoid(RP.JumpProcess, [RP.species_to_state[x] for x in S], d, solver)
end

"""
	stationary_probability_mass(MP::MarkovProcess, X::BasicSemialgebraicSet, d::Int,
							solver)

returns **lower** and **upper** bounds on the probability mass associated with the set `X`.
`d` refers to the order of the relaxation used, again the bounds will tighten
monotonically with increasing order. solver refers to the optimizer used to solve
the semidefinite programs which optimal values furnish the bounds. This is the
weakest formulation that can be used to compute bounds on the probabiltiy mass
associated with a Basic semialgebraic set. For sensible results the set `X` must have
a non-empty interior. In order to improve the bounds the user must supply a carefully
defined partition of the state space.
"""
function stationary_probability_mass(MP::MarkovProcess, X::BasicSemialgebraicSet, d::Int, solver)
	P = split_state_space(MP, X)
 	return stationary_probability_mass(MP, 1, d, solver, P)
end

function stationary_probability_mass(MP::MarkovProcess, v::Int, order::Int, solver, P::Partition)
	model, w = stationary_indicator(MP, v, order, P, solver; sense = 1) # Max
	optimize!(model)
	ub = Bound(-objective_value(model), model, P, Dict(v => value(w[v]) for v in vertices(P.graph)))
	model, w = stationary_indicator(MP, v, order, P, solver; sense = -1) # Min
	optimize!(model)
	lb = Bound(objective_value(model), model, P, Dict(v => value(w[v]) for v in vertices(P.graph)))
 	return lb, ub
end

function stationary_indicator(MP::MarkovProcess, v_target::Int, order::Int, P::Partition, solver; sense = 1)
    model = SOSModel(solver)
    @variable(model, s)
    w = Dict(v => (props(P.graph, v)[:cell] isa Singleton ?
                   @variable(model) :
                   @variable(model, [1], Poly(monomials(MP.x, 0:order)))[1]) for v in vertices(P.graph))

    for v in vertices(P.graph)
        flag = (v_target == v)
        add_stationarity_constraints!(model, MP, v, P, props(P.graph, v)[:cell], w, s + flag*sense)
    end

    for e in edges(P.graph)
        add_coupling_constraints!(model, MP, e, P, w)
    end

    @objective(model, Max, s)
    return model, w
end
