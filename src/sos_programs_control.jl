export optimal_control

"""
	optimal_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver)

returns a **lower** bound on the objective value of the (stochastic) optimal
control problem specified by CP. `μ0` encodes information about the distribution of
the initial state of the process; specifically, `μ0` maps a given monomial to the
corresponding moment of the initial distribution. trange refers to an *ordered*
set of time points discretizing the control horizon. trange[end] should coincide
with the end of the control horizon, i.e., `trange[end] = Inf` in case of an infinite
horizon problem. The bound is computed via a SOS program of degree d solved with
an appropriate method given by solver.

The bound can be tightened by populating trange or increasing `d`.
"""
function optimal_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(CP.MP.X))
	if trange[1] == 0
		trange = trange[2:end]
	end
	if isinf(trange[end])
		trange = trange[1:end-1]
	end
	if isinf(CP.T)
		model, w = infinite_horizon_control(CP, μ0, d, trange, solver, P)
	else
		model, w = finite_horizon_control(CP, μ0, d, trange, solver, P)
	end
	optimize!(model)
	return Bound(objective_value(model), model, P, dual_poly(w, CP.MP.iv, trange))
end

## Finite horizon control problems
function finite_horizon_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(CP.MP.X))
	if isa(CP.Objective, LagrangeMayer)
		model, w = finite_horizon_LM(CP, μ0, d, trange, solver, P)
	elseif isa(CP.Objective, ExitProbability)
		model, w = finite_horizon_EP(CP, μ0, d, trange, solver, P)
	elseif isa(CP.Objective, TerminalSetProbability)
		model, w = finite_horizon_TP(CP, μ0, d, trange, solver, P)
	else
		error("Objective function type not supported")
	end
	return model, w
end

## Lagrange Mayer
function finite_horizon_LM(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(CP.MP.X))
    MP = CP.MP
    t = MP.iv
    nT = length(trange)
    Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
    T = @set(t >= 0 && t <= 1)

    model = SOSModel(solver)
    w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
                    @variable(model, [1], Poly(monomials(t, 0:d)))[1] :
                    @variable(model, [1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)

    for v in vertices(P.graph)
        add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], intersect(T, CP.U), Δt, w, - CP.Objective.l)
    end

    for v in vertices(P.graph)
        add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT, v], CP.Objective.m, v)
    end

    for e in edges(P.graph)
        add_coupling_constraints!(model, MP, e, P, T, Δt, w)
    end
	obj = expectation(μ0[first(keys(μ0))] isa Dict ? [subs(w[1, v], t => 0) for v in vertices(P.graph)] : subs(w[1, 1], t => 0), μ0)

	if !isempty(CP.PathChanceConstraints)
		# TODO: improve formulation to take advantage of partition
		for PC in CP.PathChanceConstraints
			s_pc = @variable(model)
			w_pc = @variable(model, [1:nT], Poly(monomials(sort(sort(vcat(MP.x, t), rev = true),rev=true), 0:d)))
			add_path_chance_constraint!(model, MP, PC, intersect(T, CP.U), Δt, w_pc, s_pc)
			obj += s_pc*(1-PC.α)
		end
	end
	if !isempty(CP.TerminalChanceConstraints)
		for TC in CP.TerminalChanceConstraints
			s_pc = @variable(model)
			@constraint(model, s_p >= 0)
			for v in vertices(P.graph)
				add_terminal_chance_constraint!(model, MP, TC, v, P, props(P.graph, v)[:cell], w, - s_pc + CP.Objective.m)
			end
			obj += s_pc*(1-TC.α)
		end
	end
	@objective(model, Max, obj)
    return model, w
end

## Exit probability
function finite_horizon_EP(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver, P::Partition = trivial_partition(CP.MP.X))
	nT = length(trange)
	MP = CP.MP
	t = MP.iv
	Δt = [(i == 1 ? trange[i] : trange[i] - trange[i-1]) for i in 1:nT]
	T = @set(t >= 0 && t <= 1)

	model = SOSModel(solver)
	w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
                    @variable(model, [1], Poly(monomials(t, 0:d)))[1] :
                    @variable(model, [1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)

	∂X = ∂(CP.Objective.X)
	for v in vertices(P.graph)
		add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], intersect(T, CP.U), Δt, w, 0)
		add_boundary_constraints!(model, MP, v, P, props(P.graph, v)[:cell], T, w, CP.Objective.X, ∂X)
	end

	for v in vertices(P.graph)
		add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT,v], 0)
	end

	for e in edges(P.graph)
		add_coupling_constraints!(model, MP, e, P, T, Δt, w)
	end

	obj = expectation([subs(w[1,v], t => 0) for v in vertices(P.graph)], μ0)
	if !isempty(CP.PathChanceConstraints)
		# TODO: improve formulation to take advantage of partition
		for PC in CP.PathChanceConstraints
			s_pc = @variable(model)
			w_pc = @variable(model, [1:nT], Poly(monomials(sort(vcat(MP.x, t),rev=true), 0:d)))
			add_path_chance_constraint!(model, MP, PC, intersect(T, CP.U), Δt, w_pc, s_pc)
			obj += s_pc*(1-PC.α)
		end
	end
	if !isempty(CP.TerminalChanceConstraints)
		error("Terminal chance constraints currently not supported for Exit Probability problems!")
	end
	@objective(model, Max, obj)
	return model, w
end

# terminal probability
function finite_horizon_TP(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver)
	nT = length(trange)
	MP = CP.MP
	t = MP.iv
	Δt = [(i == 1 ? trange[i] : trange[i] - trange[i-1]) for i in 1:nT]
	T = @set(t >= 0 && t <= 1)

	model = SOSModel(solver)
	w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
                    @variable(model, [1], Poly(monomials(t, 0:d)))[1] :
                    @variable(model, [1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)

	for v in vertices(P.graph)
		add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], intersect(T, CP.U), Δt, w, 0)
	end

	for v in vertices(P.graph)
		add_transversality_constraints!(model, MP, intersect(props(P.graph, v)[:cell], CP.Objective.X), w[nT, v], -1)
		add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT, v], 0)
	end

	for e in edges(P.graph)
		add_coupling_constraints!(model, MP, e, P, T, Δt, w)
	end
	obj = expectation(μ0[first(keys(μ0))] isa Dict ? [subs(w[1, v], t => 0) for v in vertices(P.graph)] : subs(w[1, 1], t => 0), μ0)

	if !isempty(CP.PathChanceConstraints)
		# TODO: improve formulation to take advantage of partition
		for PC in CP.PathChanceConstraints
			s_pc = @variable(model)
			w_pc = @variable(model, [1:nT], Poly(monomials(sort(vcat(MP.x, t),rev=true), 0:d)))
			add_path_chance_constraint!(model, MP, PC, intersect(T, CP.U), Δt, w_pc, s_pc)
			obj += s_pc*(1-PC.α)
		end
	end
	if !isempty(CP.TerminalChanceConstraints)
		for TC in CP.TerminalChanceConstraints
			s_pc = @variable(model)
			@constraint(model, s_p >= 0)
			for v in vertices(P.graph)
				add_terminal_chance_constraint!(model, MP, TC, v, P, props(P.graph, v)[:cell], w, - s_pc + CP.Objective.m)
			end
			obj += s_pc*(1-TC.α)
		end
	end
	@objective(model, Max, obj)
	return model, w
end

## Discounted inf horizon control problems
function infinite_horizon_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, solver)
	@assert isa(CP.Objective, LagrangeMayer) "Objective function type not supported"
	@assert CP.Objective.m == 0 "Mayer term must be 0 for infinite horizon problem"
	@assert isempty(CP.ChancePathConstraints) "Chance path constraints are not supported in infinite horizon problems"
	@assert isempty(CP.ChanceTerminalConstraints) "Chance terminal constraints are not supported in infinite horizon problems"
	return infinite_horizon_LM(CP, μ0, d, trange, solver)
end

function infinite_horizon_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector, solver, P::Partition = trivial_partition(CP.MP.X))
    MP = CP.MP
    t = MP.iv
    nT = length(trange)
    Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
    ρ = CP.discount_factor
	T = @set(t >= 0 && t <= 1)

    model = SOSModel(solver)
	w = Dict((k, v) => (props(P.graph, v)[:cell] isa Singleton ?
                    @variable(model, [1], Poly(monomials(t, 0:d)))[1] :
                    @variable(model, [1], Poly(monomials(sort(vcat(MP.x, t), rev = true), 0:d)))[1]) for v in vertices(P.graph), k in 1:nT)
	w∞ = Dict(v => (props(P.graph, v)[:cell] isa Singleton ?
					@variable(model) :
					@variable(model, [1], Poly(monomials(MP.x, 0:d)))[1]) for v in vertices(P.graph))

	for v in vertices(P.graph)
        add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], intersect(T, CP.U), Δt, w, - CP.Objective.l; ρ = ρ)
		add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], intersect(CP.U), [1], w∞, - CP.Objective.l; ρ = ρ)
		if props(P.graph, v)[:cell] isa Singleton
			@constraint(model, w∞[v] - w[nT, v](t => 1) >= 0)
		else
			@constraint(model, w∞[v] - subs(w[nT, v], t => 1) >= 0, domain = props(P.graph, v)[:cell])
		end
		add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w∞, 0)
    end

    for e in edges(P.graph)
        add_coupling_constraints!(model, MP, e, P, T, Δt, w)
		add_coupling_constraints!(model, MP, e, P, T, [1], w∞)
    end

    @objective(model, Max, expectation([subs(w[1,v], t => 0) for v in vertices(P.graph)], μ0))
    return model, merge(w, w∞)
end

#=
function optimal_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, P::Partition, solver)
	if trange[1] == 0
		trange = trange[2:end]
	end
	if trange[end] == Inf
		trange = trange[1:end-1]
	end
	if isinf(CP.T)
        model = infinite_horizon_control(CP, μ0, d, trange, P, solver)
    else
        model = finite_horizon_control(CP, μ0, d, trange, P, solver)
    end
    optimize!(model)
    return objective_value(model), termination_status(model), MOI.get(model, MOI.SolveTime())
end

function finite_horizon_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, P::Partition, solver)
    MP = CP.MP
    t = MP.iv
    nT = length(trange)
    Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
    T = @set(t >= 0 && t <= 1)

    model = SOSModel(solver)
    w = [Dict(v => (props(P.graph, v)[:cell] isa Singleton ?
                    @variable(model, [1], Poly(monomials(t, 0:d)))[1] :
                    @variable(model, [1], Poly(monomials(vcat(MP.x, t), 0:d)))[1]) for v in vertices(P.graph))
         for k in 1:nT]

    for v in vertices(P.graph)
        add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], intersect(T, CP.U), Δt, w, - CP.Objective.l)
    end

    for v in vertices(P.graph)
        add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w[nT][v], CP.Objective.m)
    end

    for e in edges(P.graph)
        add_coupling_constraints!(model, MP, e, P, T, Δt, w)
    end

    @objective(model, Max, expectation([subs(w[1][v], t => 0) for v in vertices(P.graph)], μ0))
    return model
end

function infinite_horizon_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector, p::Partition, solver)
	if trange[1] == 0
        trange = trange[2:end]
    end
    MP = CP.MP
    t = MP.iv
    nT = length(trange)
    Δt = vcat(trange[1], [trange[i] - trange[i-1] for i in 2:nT])
    ρ = CP.discount_factor
	T = @set(t >= 0 && t <= 1)

    model = SOSModel(solver)
	w = [Dict(v => (props(P.graph, v)[:cell] isa Singleton ?
					@variable(model, [1], Poly(monomials(t, 0:d)))[1] :
					@variable(model, [1], Poly(monomials(vcat(MP.x, t), 0:d)))[1]) for v in vertices(P.graph))
		 for k in 1:nT]
	w∞ = Dict(v => (props(P.graph, v)[:cell] isa Singleton ?
					@variable(model) :
					@variable(model, [1], Poly(monomials(MP.x, 0:d)))[1]) for v in vertices(P.graph))

	for v in vertices(P.graph)
        add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], intersect(T, CP.U), Δt, w, - CP.Objective.l; ρ = ρ)
		add_dynamics_constraints!(model, MP, v, P, props(P.graph, v)[:cell], intersect(CP.U), [1], w∞, - CP.Objective.l; ρ = ρ)
		if props(P.graph, v)[:cell] isa Singleton
			@constraint(model, w∞[v] - w[nT][v](t => 1) >= 0)
		else
			@constraint(model, w∞[v] - subs(w[nT][v], t => 1) >= 0, domain = props(P.graph, v)[:cell])
		end
		add_transversality_constraints!(model, MP, props(P.graph, v)[:cell], w∞, 0)
    end

    for e in edges(P.graph)
        add_coupling_constraints!(model, MP, e, P, T, Δt, w)
		add_coupling_constraints!(model, MP, e, P, T, [1], w∞)
    end

    @objective(model, Max, expectation(w[1], μ0, P))
    return model
end
=#
