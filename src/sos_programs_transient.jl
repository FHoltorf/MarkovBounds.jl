export transient_polynomial, transient_mean, transient_variance, transient_covariance_ellipsoid

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
function transient_pop(MP::MarkovProcess, μ0::Dict, p::APL, d::Int,
                       trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(MP.X))
    if trange[1] == 0
        trange = trange[2:end]
    end
    t = MP.time
    nT = length(trange)
    Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
    T  = @set(t >= 0 && t <= 1)

    model = SOSModel(solver)
    w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
                    @variable(model, [1], Poly(monomials(t, 0:d)))[1] :
                    @variable(model, [1], Poly(monomials(vcat(MP.x, t), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)

    for v in vertices(P.graph)
        add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], T, Δt, w, 0)
    end

    for v in vertices(P.graph)
        add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT,v], p)
    end

    for e in edges(P.graph)
        add_coupling_constraints!(model, MP, e, P, T, Δt, w)
    end
	@objective(model, Max, expectation(μ0[first(keys(μ0))] isa Dict ? [subs(w[1, v], t => 0) for v in vertices(P.graph)] : subs(w[1, 1], t => 0), μ0))
    return model, w
end

function transient_pop(MP::MarkovProcess, μ0::Dict, p::Num, d::Int,
					   trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(MP.X))
	if μ0[first(keys(μ0))] isa Dict
		μ0 = Dict(v => Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[v][mono] for mono in keys(μ0[v])) for v in keys(μ0))
	else
		μ0 = Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[mono] for mono in keys(μ0))
	end
	return transient_pop(MP, μ0, polynomialize_expr(p, MP.poly_vars), d, trange, solver, P)
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
function transient_polynomial(MP::MarkovProcess, μ0::Dict, p::APL, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(MP.X))
    model, w = transient_pop(MP, μ0, p, d, trange, solver, P)
    optimize!(model)
    return Bound(objective_value(model), model, P, dual_poly(w, MP.time, trange))

end

function transient_polynomial(MP::MarkovProcess, μ0::Dict, p::Num, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(MP.X))
	if μ0[first(keys(μ0))] isa Dict
		μ0 = Dict(v => Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[v][mono] for mono in keys(μ0[v])) for v in keys(μ0))
	else
		μ0 = Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[mono] for mono in keys(μ0))
	end
	return transient_polynomial(MP, μ0, polynomialize_expr(p, MP.poly_vars), d, trange, solver, P)
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
function transient_mean(MP::MarkovProcess, μ0::Dict, p::APL, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(MP.X))
    lb = transient_polynomial(MP, μ0, p, d, trange, solver)
    ub = transient_polynomial(MP, μ0, -p, d, trange, solver)
	ub.value *= -1
    return lb, ub
end

function transient_mean(MP::MarkovProcess, μ0::Dict, p::Num, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(MP.X))
	if μ0[first(keys(μ0))] isa Dict
		μ0 = Dict(v => Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[v][mono] for mono in keys(μ0[v])) for v in keys(μ0))
	else
		μ0 = Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[mono] for mono in keys(μ0))
	end
	return transient_mean(MP, μ0, polynomialize_expr(p, MP.poly_vars), d, trange, solver, P)
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
						params::Dict = Dict(), auto_scaling = false)
	RP, S0 = reaction_process_setup(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver, params = params)
	μ0 = init_moments(RP.JumpProcess.x, S0, d + maximum(maxdegree.(RP.JumpProcess.a)))
 	return transient_mean(RP.JumpProcess, μ0, RP.species_to_state[S], d, trange, solver)
end


"""
	transient_variance(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)

returns an *upper* bound on 𝔼[v(x(T))²] - 𝔼[v(x(T))]² where v is a polynomial
and x(T) the state of the Markov process MP at time T.
"""
function transient_variance(MP::MarkovProcess, μ0::Dict, p::APL, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(MP.X))
	if trange[1] == 0
		trange = trange[2:end]
	end
	t = MP.time
	nT = length(trange)
	Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
	T  = @set(t >= 0 && t <= 1)

	model = SOSModel(solver)
	w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
                    @variable(model, [1], Poly(monomials(t, 0:d)))[1] :
                    @variable(model, [1], Poly(monomials(vcat(MP.x, t), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)
	@variable(model, S[1:2])
	@constraint(model, [1 S[1]; S[1] S[2]] in PSDCone())

	for v in vertices(P.graph)
		add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], T, Δt, w, 0)
	end

	for v in vertices(P.graph)
		add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT, v], -(p^2 + 2*S[1]*p))
	end

	for e in edges(P.graph)
		add_coupling_constraints!(model, MP, e, P, T, Δt, w)
	end

	@objective(model, Max, expectation(μ0[first(keys(μ0))] isa Dict ? [subs(w[1, v], t => 0) for v in vertices(P.graph)] : subs(w[1, 1], t => 0), μ0) - S[2])
	optimize!(model)
	return Bound(-objective_value(model), model, P, dual_poly(w, t, trange))
end

function transient_variance(MP::MarkovProcess, μ0::Dict, p::Num, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(MP.X))
	if μ0[first(keys(μ0))] isa Dict
		μ0 = Dict(v => Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[v][mono] for mono in keys(μ0[v])) for v in keys(μ0))
	else
		μ0 = Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[mono] for mono in keys(μ0))
	end
	return transient_variance(MP, μ0, polynomialize_expr(p, MP.poly_vars), d, trange, solver, P)
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
							auto_scaling = false, params::Dict = Dict())
	RP, S0 = reaction_process_setup(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver, params = params)
	μ0 = init_moments(RP.JumpProcess.x, S0, d + maximum(maxdegree.(RP.JumpProcess.a)))
 	return transient_variance(RP.JumpProcess, μ0, RP.species_to_state[S], d, trange, solver)
end

"""
	transient_covariance_ellipsoid(MP::MarkovProcess, μ0::Dict, v::Vector{APL}, d::Int, trange::AbstractVector{<:Real}, solver)

returns an *upper* bound on the volume of the covariance ellipsoid det(𝔼(v(x(T))v(x(T))ᵀ)),
where v is a polynomial and x(T) the state of the Markov process MP at time T.
"""
function transient_covariance_ellipsoid(MP::MarkovProcess, μ0::Dict, p::Vector{<:APL}, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(MP.X))
	if trange[1] == 0
		trange = trange[2:end]
	end
	t = MP.time
	nT = length(trange)
	Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
	T  = @set(t >= 0 && t <= 1)
	n = length(p)

	model = SOSModel(solver)
	w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
                    @variable(model, [1], Poly(monomials(t, 0:d)))[1] :
                    @variable(model, [1], Poly(monomials(vcat(MP.x, t), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)
	@variable(model, S[1:n+1,1:n+1], PSD)
 	@variable(model, U[1:2n,1:2n], PSD)
 	@variable(model, r[1:n])
    @variable(model, q[1:n])

	@constraint(model, [1 S[1]; S[1] S[2]] in PSDCone())
	@constraint(model, S[1:n,1:n] .== U[1:n, 1:n])
    @constraint(model, [i in 1:n, j in 1:i], -2*U[n+i,j] - (i == j ? U[n+i,n+i] - r[i] : 0) == 0)
    @constraint(model, [i in 1:n], [-1, q[i], r[i]] in MOI.DualExponentialCone())

	for v in vertices(P.graph)
		add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], T, Δt, w, 0)
	end

	for v in vertices(P.graph)
		add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT, v], -(p'*S[1:n, 1:n]*p + 2*S[end,1:n]'*p))
	end

	for e in edges(P.graph)
		add_coupling_constraints!(model, MP, e, P, T, Δt, w)
	end

	@objective(model, Max, expectation(μ0[first(keys(μ0))] isa Dict ? [subs(w[1, v], t => 0) for v in vertices(P.graph)] : subs(w[1, 1], t => 0), μ0) - S[n+1,n+1] - sum(q))
	optimize!(model)
	return Bound(exp(-objective_value(model)), model, P, dual_poly(w, t, trange))
end

function transient_covariance_ellipsoid(MP::MarkovProcess, μ0::Dict, v::Vector{Num}, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(MP.X))
	if μ0[first(keys(μ0))] isa Dict
		μ0 = Dict(v => Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[v][mono] for mono in keys(μ0[v])) for v in keys(μ0))
	else
		μ0 = Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[mono] for mono in keys(μ0))
	end
	return transient_covariance_ellipsoid(MP, μ0, polynomialize_expr(v, MP.poly_vars), d, trange, solver, P)
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
										auto_scaling = false, params::Dict = Dict())
	RP, x0 = reaction_process_setup(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver, params = params)
	μ0 = init_moments(RP.JumpProcess.x, x0, d + maximum(maxdegree.(RP.JumpProcess.a)))
 	return transient_covariance_ellipsoid(RP.JumpProcess, μ0, [RP.species_to_state[x] for x in S], d, trange, solver)
end
