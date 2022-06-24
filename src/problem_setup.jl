# for stationary problems
function add_stationarity_constraints!(model::Model, MP::MarkovProcess, v::Int, ::Partition, domain::AbstractSemialgebraicSet, w, rhs)
    @constraint(model, inf_generator(MP, w[v]) >= rhs, domain = domain)
end

function add_stationarity_constraints!(model::Model, MP::MarkovProcess, v::Int, ::Partition, domain::Vector{<:AbstractSemialgebraicSet}, w, rhs)
    cons = []
    for X in domain
        push!(cons, @constraint(model, inf_generator(MP, w[v]) >= rhs, domain = X))
    end
    return cons
end

function add_stationarity_constraints!(model::Model, MP::MarkovProcess, v::Int, P::Partition, ::Singleton, w, rhs::APL)
    @constraint(model, inf_generator(MP, w, v, P) >= rhs(MP.x => props(P.graph, v)[:cell].x))
end

function add_stationarity_constraints!(model::Model, MP::MarkovProcess, v::Int, P::Partition, ::Singleton, w, rhs)
	@constraint(model, inf_generator(MP, w, v, P) >= rhs)
end

function add_coupling_constraints!(model::Model, MP::MarkovProcess, e::Edge, P::Partition, w)
    if props(P.graph, e.dst)[:cell] isa Singleton
         @constraint(model, w[e.dst] - w[e.src](MP.x => props(P.graph, e.dst)[:cell].x) == 0)
    else
        Xs = props(P.graph, e)[:interface]
        if Xs isa AbstractVector
            for X in Xs # TODO FIX THIS
                @constraint(model, w[e.dst] - w[e.src] == 0, domain = X)
            end
        else
            @constraint(model, w[e.dst] - w[e.src] == 0, domain = Xs)
        end
    end
end

# for dynamic problems
function add_dynamics_constraints!(model::Model, MP::MarkovProcess, v::Int, P::Partition,
                                   space_domain::Singleton,
                                   time_domain::AbstractSemialgebraicSet,
                                   Δt::AbstractVector, w, rhs; ρ::Real = 0)
    nT = length(Δt)
    t = MP.iv
    model.obj_dict[Symbol("dynamics_$v")] = @constraint(model, [k in 1:nT], extended_inf_generator(MP, w, (k, v), P; scale = Δt[k]) - ρ*w[k, v] >= Δt[k]*rhs, domain = time_domain, base_name = "dynamics_$v")
    model.obj_dict[Symbol("temporal_coupling_$v")] = @constraint(model, [k in 2:nT], w[k, v](t => 0) - w[k-1, v](t => 1) >= 0, base_name = "temporal_coupling_$v")
end

function add_dynamics_constraints!(model::Model, MP::MarkovProcess, v::Int, P::Partition,
                                   space_domain::Singleton,
                                   time_domain::AbstractSemialgebraicSet,
                                   Δt::AbstractVector, w, rhs::APL; ρ::Real = 0)
    nT = length(Δt)
    t = MP.iv
    model.obj_dict[Symbol("dynamics_$v")] = @constraint(model, [k in 1:nT], extended_inf_generator(MP, w, (k, v), P; scale = Δt[k])  - ρ*w[k, v] >= Δt[k]*subs(rhs, MP.x => space_domain.x), domain = time_domain, base_name = "dynamics_$v")
    model.obj_dict[Symbol("temporal_coupling_$v")] = @constraint(model, [k in 2:nT], w[k, v](t => 0) - w[k-1,v](t => 1) >= 0, base_name = "temporal_coupling_$v")
end

function add_dynamics_constraints!(model::Model, MP::MarkovProcess, v::Int, P::Partition,
                                   space_domain::AbstractSemialgebraicSet,
                                   time_domain::AbstractSemialgebraicSet,
                                   Δt::AbstractVector, w, rhs; ρ::Real = 0)
    nT = length(Δt)
    XT = intersect(space_domain,time_domain)
    t = MP.iv
    model.obj_dict[Symbol("dynamics_$v")] = @constraint(model, [k in 1:nT], extended_inf_generator(MP, w[k, v]; scale = Δt[k]) - ρ*w[k, v] >= Δt[k]*rhs, domain = XT, base_name = "dynamics_$v")
    model.obj_dict[Symbol("temporal_coupling_$v")] = @constraint(model, [k in 2:nT], subs(w[k, v], t => 0) - subs(w[k-1, v], t => 1) >= 0, domain = space_domain, base_name = "temporal_coupling_$v")
end

function add_dynamics_constraints!(model::Model, MP::MarkovProcess, v::Int, P::Partition,
                                   space_domain::Vector{<:AbstractSemialgebraicSet},
                                   time_domain::AbstractSemialgebraicSet,
                                   Δt::AbstractVector, w, rhs; ρ::Real = 0)
    nT = length(Δt)
    t = MP.iv
    model.obj_dict[Symbol("dynamics_$v")] = @constraint(model, [k in 1:nT, X in space_domain], extended_inf_generator(MP, w[k, v]; scale = Δt[k]) - ρ*w[k, v] >= Δt[k]*rhs, domain = XT = intersect(X,time_domain), base_name = "dynamics_$v")
    model.obj_dict[Symbol("temporal_coupling_$v")] = @constraint(model, [k in 2:nT, X in space_domain], subs(w[k, v], t => 0) - subs(w[k-1,v], t => 1) >= 0, domain = X, base_name = "temporal_coupling_$v")
end

# for infinite horizon problems
function add_tail_constraints!(model::Model, MP::MarkovProcess, v::Int, P::Partition,
                               space_domain::Singleton,
                               control_domain::AbstractSemialgebraicSet,
                               w∞, rhs; ρ::Real = 0)
    model.obj_dict[Symbol("tail_dynamics_$v")] = @constraint(model, inf_generator(MP, w∞, v, P) - ρ*w∞[v] >= rhs, domain = control_domain, base_name = "tail_dynamics_$v")
end

function add_tail_constraints!(model::Model, MP::MarkovProcess, v::Int, P::Partition,
                               space_domain::Singleton,
                               control_domain::AbstractSemialgebraicSet,
                               w∞, rhs::APL; ρ::Real = 0)
    model.obj_dict[Symbol("tail_dynamics_$v")] = @constraint(model, inf_generator(MP, w∞, v, P)  - ρ*w∞[v] >= rhs(MP.x => space_domain.x), domain = control_domain, base_name = "tail_dynamics_$v")
end

function add_tail_constraints!(model::Model, MP::MarkovProcess, v::Int, P::Partition,
                               space_domain::AbstractSemialgebraicSet,
                               control_domain::AbstractSemialgebraicSet,
                               w∞, rhs; ρ::Real = 0)
    model.obj_dict[Symbol("tail_dynamics_$v")] = @constraint(model, inf_generator(MP, w∞[v]) - ρ*w∞[v] >= rhs, domain = intersect(space_domain, control_domain), base_name = "tail_dynamics_$v")
end

function add_tail_constraints!(model::Model, MP::MarkovProcess, v::Int, P::Partition,
                               space_domain::Vector{<:AbstractSemialgebraicSet},
                               control_domain::AbstractSemialgebraicSet,
                               w∞, rhs; ρ::Real = 0)
    model.obj_dict[Symbol("tail_dynamics_$v")] = @constraint(model, [X in space_domain], inf_generator(MP, w∞[v]) - ρ*w∞[v] >= rhs, domain = intersect(X, control_domain), base_name = "tail_dynamics_$v")
end

# terminal cost transversality constraints
function add_transversality_constraints!(model::Model, MP::MarkovProcess, space_domain::Singleton, w, rhs::Real, v::Int)
    model.obj_dict[Symbol("transversality_$v")] = @constraint(model, w(MP.iv => 1) <= rhs, base_name = "transversality_$v")
end

function add_transversality_constraints!(model::Model, MP::MarkovProcess, space_domain::Singleton, w, rhs::APL, v::Int)
    model.obj_dict[Symbol("transversality_$v")] = @constraint(model, w(MP.iv => 1) <= rhs(MP.x => space_domain.x), base_name = "transversality_$v")
end

function add_transversality_constraints!(model::Model, MP::MarkovProcess, space_domain::AbstractSemialgebraicSet, w, rhs, v::Int)
    model.obj_dict[Symbol("transversality_$v")] = @constraint(model, subs(w, MP.iv => 1) <= rhs, domain = space_domain, base_name = "transversality_$v")
end

function add_transversality_constraints!(model::Model, MP::MarkovProcess, space_domain::Singleton, w, rhs, v::Int)
    model.obj_dict[Symbol("transversality_$v")] = @constraint(model, w(MP.iv => 1) <= rhs, base_name = "transversality_$v")
end

function add_transversality_constraints!(model::Model, MP::MarkovProcess, space_domain::Vector{<:AbstractSemialgebraicSet}, w, rhs, v::Int)
    model.obj_dict[Symbol("transversality_$v")] = @constraint(model, [X in space_domain], subs(w, MP.iv => 1) <= rhs, domain = X, base_name = "transversality_$v")
end

function add_coupling_constraints!(model::Model, MP::MarkovProcess, e::Edge, P::Partition,
                                   time_domain::AbstractSemialgebraicSet, Δt::AbstractVector, w)
    nT = length(Δt)
    if props(P.graph, e.dst)[:cell] isa Singleton
        @constraint(model, [k in 1:nT], w[k, e.dst] - subs(w[k, e.src], MP.x => props(P.graph, e.dst)[:cell].x) == 0, domain = time_domain)
    else
        Xs = props(P.graph, e)[:interface]
        for X in Xs
            XT = intersect(X,time_domain)
            @constraint(model, [k in 1:nT], w[k, e.dst] - w[k, e.src] == 0, domain = XT)
        end
    end
end

function add_path_chance_constraint!(model::Model, MP::MarkovProcess, PC::ChanceConstraint, time_domain, Δt::AbstractVector, w_pc, s_pc)
	∂X = ∂(PC.X)
	t = MP.iv
	nT = length(Δt)
	@constraint(model, [i in 1:nT], extended_inf_generator(MP, w_pc[i]; scale = Δt[i]) >= 0, domain = intersect(PC.X, time_domain))
	@constraint(model, [i in 2:nT], subs(w_pc[i], t => 0) - subs(w_pc[i-1], t => 1) >= 0, domain = PC.X)
	@constraint(model, [i in 1:nT, k in 1:length(∂X)], w_pc[i] >= 0, domain = ∂X[k])
	@constraint(model, -s_pc - subs(w_pc[nT], t => 1) >= 0, domain = PC.X)
	@constraint(model, s_pc >= 0)
end

function add_boundary_constraints!(model::Model, MP::MarkovProcess, nT::Int, v::Int, P::Partition, space_domain::AbstractSemialgebraicSet,
								   time_domain::AbstractSemialgebraicSet, w, X::AbstractSemialgebraicSet, ∂X::AbstractVector{<:AbstractSemialgebraicSet})
	@constraint(model, [i in 1:nT, k in 1:length(∂X)], w[i, v] >= 0, domain = intersect(∂X[k], space_domain, T))
	@constraint(model, - 1 - subs(w[nT, v], MP.iv => 1) >= 0, domain = intersect(space_domain, X))
end

function add_terminal_chance_constraint!(model::Model, MP::MarkovProcess, TC::ChanceConstraint, v::Int, P::Partition, space_domain::AbstractSemialgebraicSet, w, rhs)
	@constraint(model, subs(w[nT, v], MP.iv => 1) <= rhs, domain = intersect(TC.X,space_domain))
end
