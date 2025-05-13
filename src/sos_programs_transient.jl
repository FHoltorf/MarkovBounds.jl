export transient_polynomial, transient_mean, transient_variance, transient_covariance_ellipsoid

@doc raw"""
	transient_pop(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)

returns SOS program of degree `d` for computing a **lower** bound on ``\mathbb{E}[v(x(T))]``
where ``v`` is a polynomial and ``x(T)`` the state of the Markov process `MP` at time `T = trange[end]`.
μ0 encodes the distribution of the initial state of the process in terms of its moments;
specifically, it maps monomials to the respective moments of the initial distribution.
`trange` is an **ordered** collection of time points used to discretize the time
horizon ``[0,T]``, i.e., `T = trange[end]`. Populating `trange` and increasing `d` has a
tightening effect on the bound furnished by the assembled SOS program.
"""
function transient_pop(MP::MarkovProcess, μ0::Dict, p::APL, d::Int,
                       trange::AbstractVector{<:Real}, solver, 
					   P::Partition = trivial_partition(MP.X);
					   inner_approx = SOSCone)
    if trange[1] == 0
        trange = trange[2:end]
    end
    t = MP.iv
    nT = length(trange)
    Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
    T  = @set(t >= 0 && t <= 1)

    model = SOSModel(solver)
	PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, inner_approx)
    w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
                    @variable(model, [1:1], Poly(monomials(t, 0:d)))[1] :
                    @variable(model, [1:1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)

    for v in vertices(P.graph)
        add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], T, Δt, w, 0)
    end

    for v in vertices(P.graph)
        add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT,v], p, v)
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

@doc raw"""
	transient_polynomial(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)

returns a **lower** bound on ``\mathbb{E}[v(x(T))]`` where ``v`` is a polynomial and ``x(T)`` the state
of the Markov process `MP` at time `T = trange[end]`. `μ0` encodes the
distribution of the initial state of the process in terms of its moments; specifically,
it maps monomials to the respective moments of the initial distribution.
`trange` is an **ordered** collection of time points used to discretize the time
horizon ``[0,T]``, i.e., `T = trange[end]`. Populating `trange` and increasing `d` improves
the computed bound.
"""
function transient_polynomial(MP::MarkovProcess, μ0::Dict, p::APL, d::Int, trange::AbstractVector{<:Real}, solver, 
							 P::Partition = trivial_partition(MP.X);
							 inner_approx = SOSCone)
    model, w = transient_pop(MP, μ0, p, d, trange, solver, P, inner_approx=inner_approx)
    optimize!(model)
    return Bound(objective_value(model), model, P, dual_poly(w, MP.iv, trange))
end

function transient_polynomial(MP::MarkovProcess, μ0::Dict, p::Num, d::Int, trange::AbstractVector{<:Real}, solver, 
							  P::Partition = trivial_partition(MP.X);
							  inner_approx = SOSCone)
	if μ0[first(keys(μ0))] isa Dict
		μ0 = Dict(v => Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[v][mono] for mono in keys(μ0[v])) for v in keys(μ0))
	else
		μ0 = Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[mono] for mono in keys(μ0))
	end
	return transient_polynomial(MP, μ0, polynomialize_expr(p, MP.poly_vars), d, trange, solver, P, inner_approx=inner_approx)
end

@doc raw"""
	transient_mean(MP::MarkovProcess, μ0::Dict, x::APL, d::Int, trange::AbstractVector{<:Real}, solver)

returns a **lower** and **upper** bound on ``\mathbb{E}[v(x(T))]`` where ``v`` is a polynomial and ``x(T)`` the state
of the Markov process `MP` at time `T = trange[end]`. `μ0` encodes the
distribution of the initial state of the process in terms of its moments; specifically,
it maps monomials to the respective moments of the initial distribution.
trange is an **ordered** collection of time points used to discretize the time
horizon ``[0,T]``, i.e., `T = trange[end]`. Populating `trange` and increasing `d` improves
the computed bounds.
"""
function transient_mean(MP::MarkovProcess, μ0::Dict, p::APL, d::Int, trange::AbstractVector{<:Real}, solver, 
						P::Partition = trivial_partition(MP.X);
						inner_approx = SOSCone)
    lb = transient_polynomial(MP, μ0, p, d, trange, solver, P, inner_approx=inner_approx)
    ub = transient_polynomial(MP, μ0, -p, d, trange, solver, P, inner_approx=inner_approx)
	ub.value *= -1
    return lb, ub
end

function transient_mean(MP::MarkovProcess, μ0::Dict, p::Num, d::Int, trange::AbstractVector{<:Real}, solver, 
						P::Partition = trivial_partition(MP.X);
						inner_approx = SOSCone)
	if μ0[first(keys(μ0))] isa Dict
		μ0 = Dict(v => Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[v][mono] for mono in keys(μ0[v])) for v in keys(μ0))
	else
		μ0 = Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[mono] for mono in keys(μ0))
	end
	return transient_mean(MP, μ0, polynomialize_expr(p, MP.poly_vars), d, trange, solver, P, inner_approx=inner_approx)
end

"""
	transient_mean(rn::ReactionSystem, S0::Dict, S, d::Int, trange::AbstractVector{<:Number}, solver,
			scales = Dict(s => 1 for s in species(rn));
			auto_scaling = false)

returns a **lower** and **upper** bound on the mean of the molecular count of species
`S` in reaction network `rn` at time `T = trange[end]`. `S0` refers to the **deterministic** initial state
of the reaction system (including all species!).
`trange` is an **ordered** collection of time points used to discretize the time
horizon ``[0,T]``, i.e., `T = trange[end]`. Populating trange and increasing `d` improves
the computed bounds.

For numerical stability, it is recommended to provide scales of the expected
magnitude of molecular counts for the different species at steady state. If
the system is **closed** it is also possible to enable auto_scaling which will
find the maximum molecular counts for each species under stoichiometry
constraints (via LP).
"""
function transient_mean(rn::ReactionSystem, S0::Dict, S, d::Int, trange::AbstractVector{<:Number}, solver,
		  				scales = Dict(s => 1 for s in species(rn));
						params::Dict = Dict(), auto_scaling = false, 
						inner_approx = SOSCone)
	RP, S0 = reaction_process_setup(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver, params = params)
	μ0 = init_moments(RP.JumpProcess.x, S0, d + maximum(maxdegree.(RP.JumpProcess.a)))
 	return transient_mean(RP.JumpProcess, μ0, RP.species_to_state[S], d, trange, solver, inner_approx=inner_approx)
end


@doc raw"""
	transient_variance(MP::MarkovProcess, μ0::Dict, v::APL, d::Int, trange::AbstractVector{<:Real}, solver)

returns an **upper** bound on ``\mathbb{E}[v(x(T))^2] - \mathbb{E}[v(x(T))]^2`` where ``v`` is a polynomial
and ``x(T)`` the state of the Markov process `MP` at time `T = trange[end]`.
"""
function transient_variance(MP::MarkovProcess, μ0::Dict, p::APL, d::Int, trange::AbstractVector{<:Real}, solver, 
							P::Partition = trivial_partition(MP.X);
							inner_approx = SOSCone)
	if trange[1] == 0
		trange = trange[2:end]
	end
	t = MP.iv
	nT = length(trange)
	Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
	T  = @set(t >= 0 && t <= 1)
	
	model = SOSModel(solver)
	PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, inner_approx)
	w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
                    @variable(model, [1:1], Poly(monomials(t, 0:d)))[1] :
                    @variable(model, [1:1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)
	@variable(model, S[1:2])
	@constraint(model, [1 S[1]; S[1] S[2]] in PSDCone())

	for v in vertices(P.graph)
		add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], T, Δt, w, 0)
	end

	for v in vertices(P.graph)
		add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT, v], -(p^2 + 2*S[1]*p), v)
	end

	for e in edges(P.graph)
		add_coupling_constraints!(model, MP, e, P, T, Δt, w)
	end

	@objective(model, Max, expectation(μ0[first(keys(μ0))] isa Dict ? [subs(w[1, v], t => 0) for v in vertices(P.graph)] : subs(w[1, 1], t => 0), μ0) - S[2])
	optimize!(model)
	return Bound(-objective_value(model), model, P, dual_poly(w, t, trange))
end

function transient_variance(MP::MarkovProcess, μ0::Dict, p::Num, d::Int, trange::AbstractVector{<:Real}, solver, 
							P::Partition = trivial_partition(MP.X);
							inner_approx = SOSCone)
	if μ0[first(keys(μ0))] isa Dict
		μ0 = Dict(v => Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[v][mono] for mono in keys(μ0[v])) for v in keys(μ0))
	else
		μ0 = Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[mono] for mono in keys(μ0))
	end
	return transient_variance(MP, μ0, polynomialize_expr(p, MP.poly_vars), d, trange, solver, P, inner_approx=inner_approx)
end

"""
	transient_variance(rn::ReactionSystem, S0::Dict, S, d::Int, trange::AbstractVector{<:Real}, solver,
			    scales = Dict(s => 1 for s in species(rn));
				auto_scaling = false)

returns an **upper** bound on the variance of species `S` in the reaction network `rn`
at time `T = trange[end]`. `S0` refers to the **deterministic** initial state
of the reaction system (including all species!).
`trange` is an **ordered** collection of time points used to discretize the time
horizon ``[0,T]``, i.e., `T = trange[end]`. Populating `trange` and increasing `d` improves
the computed bounds.

For numerical stability, it is recommended to provide scales of the expected
magnitude of molecular counts for the different species at steady state. If
the system is **closed** it is also possible to enable `auto_scaling` which will
find the maximum molecular counts for each species under stoichiometry
constraints (via LP).
"""
function transient_variance(rn::ReactionSystem, S0::Dict, S, d::Int, trange::AbstractVector{<:Real}, solver,
						    scales = Dict(s => 1 for s in species(rn));
							auto_scaling = false, params::Dict = Dict(),
							inner_approx = SOSCone)
	RP, S0 = reaction_process_setup(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver, params = params)
	μ0 = init_moments(RP.JumpProcess.x, S0, d + maximum(maxdegree.(RP.JumpProcess.a)))
 	return transient_variance(RP.JumpProcess, μ0, RP.species_to_state[S], d, trange, solver, inner_approx=inner_approx)
end

@doc raw"""
	transient_covariance_ellipsoid(MP::MarkovProcess, μ0::Dict, v::Vector{APL}, d::Int, trange::AbstractVector{<:Real}, solver)

returns an **upper** bound on the volume of the covariance ellipsoid ``\text{det}(\mathbb{E}[v(x(T))v(x(T))^\top] - \mathbb{E}[v(x(T))] \mathbb{E}[v(x(T))]^\top)``,
where ``v`` is a polynomial and ``x(T)`` the state of the Markov process `MP` at time `T = trange[end]`.
"""
function transient_covariance_ellipsoid(MP::MarkovProcess, μ0::Dict, p::Vector{<:APL}, d::Int, trange::AbstractVector{<:Real}, solver, 
										P::Partition = trivial_partition(MP.X);
										inner_approx = SOSCone)
	if trange[1] == 0
		trange = trange[2:end]
	end
	t = MP.iv
	nT = length(trange)
	Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
	T  = @set(t >= 0 && t <= 1)
	n = length(p)

	model = SOSModel(solver)
	PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, inner_approx)
	w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
                    @variable(model, [1:1], Poly(monomials(t, 0:d)))[1] :
                    @variable(model, [1:1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)
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
		add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT, v], -(p'*S[1:n, 1:n]*p + 2*S[end,1:n]'*p), v)
	end

	for e in edges(P.graph)
		add_coupling_constraints!(model, MP, e, P, T, Δt, w)
	end

	@objective(model, Max, expectation(μ0[first(keys(μ0))] isa Dict ? [subs(w[1, v], t => 0) for v in vertices(P.graph)] : subs(w[1, 1], t => 0), μ0) - S[n+1,n+1] - sum(q))
	optimize!(model)
	return Bound(exp(-objective_value(model)), model, P, dual_poly(w, t, trange))
end

function transient_covariance_ellipsoid(MP::MarkovProcess, μ0::Dict, v::Vector{Num}, d::Int, trange::AbstractVector{<:Real}, solver,
										P::Partition = trivial_partition(MP.X);
										inner_approx = SOSCone)
	if μ0[first(keys(μ0))] isa Dict
		μ0 = Dict(v => Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[v][mono] for mono in keys(μ0[v])) for v in keys(μ0))
	else
		μ0 = Dict(polynomialize_expr(mono, MP.poly_vars) => μ0[mono] for mono in keys(μ0))
	end
	return transient_covariance_ellipsoid(MP, μ0, polynomialize_expr(v, MP.poly_vars), d, trange, solver, P, inner_approx=inner_approx)
end

"""
	transient_covariance_ellipsoid(rn::ReactionSystem, S0::Dict, S::AbstractVector, d::Int, trange::AbstractVector{<:Real}, solver,
							scales = Dict(s => 1 for s in species(rn));
							auto_scaling = false)

returns an **upper** bound on the volume of the covariance ellipsoid associated with any
collection of chemical species in the reaction network `rn` at time `T = trange[end]`.
S0 refers to the **deterministic** initial state of the reaction system
(including all species!). `trange` is an **ordered** collection of time points used to discretize the time
horizon ``[0,T]``, i.e., `T = trange[end]`. Populating `trange` and increasing `d` improves
the computed bounds.

For numerical stability, it is recommended to provide scales of the expected
magnitude of molecular counts for the different species at steady state. If
the system is **closed** it is also possible to enable auto_scaling which will
find the maximum molecular counts for each species under stoichiometry
constraints (via LP).
"""
function transient_covariance_ellipsoid(rn::ReactionSystem, S0::Dict, S::AbstractVector, d::Int, trange::AbstractVector{<:Real}, solver,
										scales = Dict(s => 1 for s in species(rn));
										auto_scaling = false, params::Dict = Dict(),
										inner_approx = SOSCone)
	RP, x0 = reaction_process_setup(rn, S0, scales = scales, auto_scaling = auto_scaling, solver = solver, params = params)
	μ0 = init_moments(RP.JumpProcess.x, x0, d + maximum(maxdegree.(RP.JumpProcess.a)))
 	return transient_covariance_ellipsoid(RP.JumpProcess, μ0, [RP.species_to_state[x] for x in S], d, trange, solver, inner_approx=inner_approx)
end


# distributed 
function transient_indicator(MP::MarkovProcess, μ0::Dict, v_target::Int, d::Int,
							 trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(MP.X); 
							 sense = 1, inner_approx = SOSCone)
	return transient_indicator(MP, μ0, [v_target], d, trange, solver, P, sense = sense, inner_approx=inner_approx)
end

function transient_indicator(MP::MarkovProcess, μ0::Dict, v_target::AbstractArray{Int}, d::Int,
							trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(MP.X); 
							sense = 1, inner_approx = SOSCone)
	if trange[1] == 0
		trange = trange[2:end]
	end
	t = MP.iv
	nT = length(trange)
	Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
	T  = @set(t >= 0 && t <= 1)

	model = SOSModel(solver)
	PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, inner_approx)
	w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
						@variable(model, [1:1], Poly(monomials(t, 0:d)))[1] :
						@variable(model, [1:1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)

	for v in vertices(P.graph)
		add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], T, Δt, w, 0)
	end

	for v in vertices(P.graph)
		flag = (v in v_target)
		add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT,v], sense*flag, v)
	end

	for e in edges(P.graph)
		add_coupling_constraints!(model, MP, e, P, T, Δt, w)
	end
	@objective(model, Max, expectation(μ0[first(keys(μ0))] isa Dict ? [subs(w[1, v], t => 0) for v in vertices(P.graph)] : subs(w[1, 1], t => 0), μ0))
	
	optimize!(model)
	return Bound(objective_value(model), model, P, dual_poly(w, t, trange))
end

# add transient probability mass
function approximate_transient_measure(MP::MarkovProcess, μ0::Dict, p::APL, d::Int,
									   trange::AbstractVector{<:Real}, solver, 
									   P::Partition = trivial_partition(MP.X);
									   inner_approx = SOSCone)
	if trange[1] == 0
		trange = trange[2:end]
	end
	t = MP.iv
	nT = length(trange)
	Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
	T  = @set(t >= 0 && t <= 1)

	model = SOSModel(solver)
	PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, inner_approx)
	w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
						@variable(model, [1:1], Poly(monomials(t, 0:d)))[1] :
						@variable(model, [1:1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)

	for v in vertices(P.graph)
		add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], T, Δt, w, 0)
	end

	cons = Dict()
	for v in vertices(P.graph)
		cons[v] = add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT,v], p, v)
	end

	for e in edges(P.graph)
		add_coupling_constraints!(model, MP, e, P, T, Δt, w)
	end
	@objective(model, Max, expectation(μ0[first(keys(μ0))] isa Dict ? [subs(w[1, v], t => 0) for v in vertices(P.graph)] : subs(w[1, 1], t => 0), μ0))
	
	optimize!(model)
	dist = []
    for v in vertices(P.graph)
        if props(P.graph, v)[:cell] isa Singleton
            push!(dist, dual(cons[v]).a[end])
        elseif props(P.graph, v)[:cell] isa Vector{BasicSemialgebraicSet}
            push!(dist, sum(dual.(cons[v])).a[end])
        elseif props(P.graph, v)[:cell] isa BasicSemialgebraicSet
            push!(dist, dual(cons[v]).a[end])
        end
    end
	return Bound(objective_value(model), model, P, dual_poly(w, t, trange)), dist
end
"""
    max_entropy_measure(MP::MarkovProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition;
                        side_infos::BasicSemialgebraicSet)

returns approximate values for the transient measure which maximizes the entropy on the partition `P` at time t.  
"""
function max_entropy_measure(MP::MarkovProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition; inner_approx = SOSCone)
	if trange[1] == 0
		trange = trange[2:end]
	end
	t = MP.iv
	nT = length(trange)
	Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
	T  = @set(t >= 0 && t <= 1)

	model = SOSModel(solver)
	PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, inner_approx)
	w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
						@variable(model, [1:1], Poly(monomials(t, 0:d)))[1] :
						@variable(model, [1:1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)

	for v in vertices(P.graph)
		add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], T, Δt, w, 0)
	end

	cons = Dict()
	@variable(model, u[vertices(P.graph)])
	@variable(model, q[vertices(P.graph)])
	for v in vertices(P.graph)
		cons[v] = add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT,v], -q[v], v)
	end

	for e in edges(P.graph)
		add_coupling_constraints!(model, MP, e, P, T, Δt, w)
	end

	@constraint(model, [v in vertices(P.graph)], [-1, q[v], u[v]] in MOI.DualExponentialCone())

	@objective(model, Max, expectation(μ0[first(keys(μ0))] isa Dict ? [subs(w[1, v], t => 0) for v in vertices(P.graph)] : subs(w[1, 1], t => 0), μ0) 
						   - sum(u))
	
	optimize!(model)
	dist = []
    for v in vertices(P.graph)
        if props(P.graph, v)[:cell] isa Singleton
            push!(dist, dual(cons[v]).a[end])
        elseif props(P.graph, v)[:cell] isa Vector{BasicSemialgebraicSet}
            push!(dist, sum(dual.(cons[v])).a[end])
        elseif props(P.graph, v)[:cell] isa BasicSemialgebraicSet
            push!(dist, dual(cons[v]).a[end])
        end
    end
	return Bound(objective_value(model), model, P, dual_poly(w, t, trange)), dist
end

function exit_pop(MP::MarkovProcess, μ0::Dict, p::APL, ∂P, d::Int,
					trange::AbstractVector{<:Real}, solver,
					P::Partition = trivial_partition(MP.X);
					inner_approx = SOSCone)
	return exit_pop(MP, μ0, p, p, ∂P, d,
					trange, solver, P;
					inner_approx = inner_approx)
end
function exit_pop(MP::MarkovProcess, μ0::Dict, p::APL, q::APL, ∂P, d::Int,
                  trange::AbstractVector{<:Real}, solver,
				  P::Partition = trivial_partition(MP.X);
				  inner_approx = SOSCone)
    if trange[1] == 0
        trange = trange[2:end]
    end
    t = MP.iv
    nT = length(trange)
    Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
    T  = @set(t >= 0 && t <= 1)
	if length(trange) > 1
		ps = vcat(subs(p, MP.iv => Δt[1]*MP.iv), [subs(p, MP.iv => Δt[i]*MP.iv + trange[i-1]) for i in 2:nT])
	else
		ps = [subs(p, MP.iv => Δt[1]*MP.iv)]
	end
    model = SOSModel(solver)
	PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, inner_approx)
    w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
		@variable(model, [1:1], Poly(monomials(t, 0:d)))[1] :
		@variable(model, [1:1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)

    for v in vertices(P.graph)
        add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], T, Δt, w, 0)
    end
	for e in edges(P.graph)
        add_coupling_constraints!(model, MP, e, P, T, Δt, w)
    end

	# spatial exists
    for v in keys(∂P)
		if !isempty(∂P[v])
        	add_exit_constraints!(model, nT, v, props(P.graph, v)[:cell], ∂P[v], T, w, ps)
		end
    end
	

	# temporal exits
	for v in vertices(P.graph)
		add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT,v], subs(q, MP.iv => trange[end]), v)
	end

	@objective(model, Max, expectation(μ0[first(keys(μ0))] isa Dict ? [subs(w[1, v], t => 0) for v in vertices(P.graph)] : subs(w[1, 1], t => 0), μ0))

	optimize!(model)
    return Bound(objective_value(model), model, P, dual_poly(w, t, trange))
end

function exit_pop_unbounded(MP::MarkovProcess, μ0::Dict, p::APL, ∂P, d::Int,
                  trange::AbstractVector{<:Real}, solver,
				  P::Partition = trivial_partition(MP.X);
				  inner_approx = SOSCone)
    if trange[1] == 0
        trange = trange[2:end]
    end
    t = MP.iv
    nT = length(trange)
    Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
    T  = @set(t >= 0 && t <= 1)
	T_final = @set(t >= 0)
	if length(trange) > 1
		ps = vcat(subs(p, MP.iv => Δt[1]*MP.iv), [subs(p, MP.iv => Δt[i]*MP.iv + trange[i-1]) for i in 2:nT])
	else
		ps = [subs(p, MP.iv => Δt[1]*MP.iv)]
	end
    model = SOSModel(solver)
	PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, inner_approx)
    w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
		@variable(model, [1:1], Poly(monomials(t, 0:d)))[1] :
		@variable(model, [1:1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)
	w_final = Dict(v => (props(P.graph, v)[:cell] isa Singleton ?
					@variable(model, [1:1], Poly(monomials(t, 0:d)))[1] :
					@variable(model, [1:1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph))
    for v in vertices(P.graph)
        add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], T, Δt, w, 0)

		# last element 
		XT = intersect(props(P.graph, v)[:cell],T_final)
		t = MP.iv
		@constraint(model, [1:1], extended_inf_generator(MP, w_final[v]) >= 0.0, domain = XT)
		@constraint(model, [1:1], subs(w_final[v], t => 0) - subs(w[length(Δt), v], t => 1) >= 0, domain = props(P.graph, v)[:cell])
    end
	for e in edges(P.graph)
        add_coupling_constraints!(model, MP, e, P, T, Δt, w)

		# last element
		if props(P.graph, e.dst)[:cell] isa Singleton
			@constraint(model, [1:1], w_final[e.dst] - subs(w_final[e.src], MP.x => props(P.graph, e.dst)[:cell].x) == 0, domain = T_final)
		else
			Xs = props(P.graph, e)[:interface]
			for X in Xs
				XT = intersect(X,T_final)
				@constraint(model, [1:1], w_final[e.dst] - w_final[e.src] == 0, domain = XT)
			end
		end
    end

	# spatial exists
    for v in keys(∂P)
		if !isempty(∂P[v])
        	add_exit_constraints!(model, nT, v, props(P.graph, v)[:cell], ∂P[v], T, w, ps)

			# last element
			X = props(P.graph, v)[:cell]
			∂X = ∂P[v]
			ps_final = subs(p, MP.iv => trange[end] + MP.iv)
			@constraint(model, [k in 1:length(∂X)], w_final[k] <= ps_final, domain = intersect(∂X[k], X, T_final))
		end
    end
	@objective(model, Max, expectation(μ0[first(keys(μ0))] isa Dict ? [subs(w[1, v], t => 0) for v in vertices(P.graph)] : subs(w[1, 1], t => 0), μ0))

	optimize!(model)
    return Bound(objective_value(model), model, P, dual_poly(w, t, trange))
end

