using LightGraphs: edges, vertices, SimpleGraph, add_edge!, Edge
using MetaGraphs: MetaGraph, props, set_prop!

export grid_graph, Partition, props

struct Partition
    graph
    get_vertex
end

function grid_graph(x, lb, ub, n; inf_top = zeros(Int64, length(ub)), inf_floor = zeros(Int64, length(lb)))
    @assert all(lb .<= ub)

    x_ranges = []
    for k in 1:length(n)
        n_eff = n[k] - inf_top[k] - inf_floor[k]
        if n_eff < 0
            push!(x_ranges, [inf_floor[k] == 0 ? lb[k] : -Inf, (inf_top[k] == 0 ? ub[k] : Inf)])
        elseif n_eff == 0
            push!(x_ranges, [inf_floor[k] == 0 ? lb[k] : -Inf, (lb[k] + ub[k])/2, (inf_top[k] == 0 ? ub[k] : Inf)])
        else
            push!(x_ranges, [inf_floor[k] == 0 ? lb[k] : -Inf, range(lb[k], stop=ub[k], length=n_eff)..., (inf_top[k] == 0 ? ub[k] : Inf)])
        end
    end
    return grid_graph(x, x_ranges)
end

function grid_graph(x, x_ranges)
    n = length.(x_ranges) .- 1
    nv = prod(n)

    mg = MetaGraph(nv)
    for i in 1:nv
        idx = invert_index(i, n)
        subset = FullSpace()
        for k in 1:length(idx)
            if !isinf(x_ranges[k][idx[k]])
                subset = intersect(subset, @set(x[k] >= x_ranges[k][idx[k]]))
            end
            if !isinf(x_ranges[k][idx[k]+1])
                subset = intersect(subset, @set(x[k] <= x_ranges[k][idx[k]+1]))
            end
        end
        set_prop!(mg, i, :cell, subset)
        for k in 1:length(idx)
            if idx[k] < n[k]
                idx[k] += 1
                j = linearize_index(idx, n)
                add_edge!(mg, i, j)
                set_prop!(mg, Edge(i,j), :interface, (k, x_ranges[k][idx[k]]))
                idx[k] -= 1
            end
        end
    end
    get_vertex = function (x)
                    idx = zeros(Int64, length(n))
                    for k in 1:length(n)
                        j = findfirst(m -> m >= x[k], x_ranges[k]) - 1
                        if isnothing(j)
                            error("state outside domain")
                        end
                        idx[k] = max(j, 1)
                    end
                    return linearize_index(idx, n)
                 end
    return mg, get_vertex
end

function optimal_control(CP::ControlProcess, μ0::Dict, d::Int, trange::AbstractVector{<:Real}, p::Partition, solver)
    if isinf(CP.T)
        model = infinite_horizon_control(CP, μ0, d, trange, p, solver)
    else
        model = finite_horizon_control(CP, μ0, d, trange, p, solver)
    end
    optimize!(model)
    return objective_value(model), termination_status(model), MOI.get(model, MOI.SolveTime())
end

function finite_horizon_control(CP::ControlProcess, μ0::Dict, order::Int, trange::AbstractVector{<:Real}, p::Partition, solver)
    MP = CP.MP
    @assert typeof(MP) == DiffusionProcess
    t = MP.time
    nₜ = length(trange)
    Δt = [trange[1], [trange[i] - trange[i-1] for i in 2:nₜ]...]

    model = SOSModel(solver)

    @variable(model, w[vertices(p.graph), k in 1:nₜ], Poly(monomials([MP.x..., t], 0:order)))

    for v in vertices(p.graph), k in 1:nₜ
        X = intersect(props(p.graph, v)[:cell], @set(t >= 0 && 1-t >= 0), CP.U)
        @constraint(model, extended_inf_generator(MP, w[v,k]; scale = Δt[k]) + Δt[k]*CP.Objective.l >= 0, domain = X)
    end

    for v in vertices(p.graph), k in 2:nₜ
        @constraint(model, subs(w[v,k], t => 0) - subs(w[v,k-1], t => 1) >= 0, domain=props(p.graph, v)[:cell])
    end

    for e in edges(p.graph), k in 1:nₜ
        i, val = props(p.graph, e)[:interface]
        if all(subs(MP.σ, MP.x[i] => val) .== 0)
            X = FullSpace()
            for p in props(p.graph, e.dst)[:cell].p
                if !(variables(p) == [MP.x[i]])
                    p = subs(p, MP.x[i] => val)
                    X = intersect(X, @set(p >= 0))
                end
            end
            X = intersect(X, CP.U, @set(t >= 0 && 1-t >= 0))
            @constraint(model, subs((w[e.dst,k] - w[e.src,k])*MP.f[i], MP.x[i] => val) >= 0, domain=X)
        else
            @constraint(model, subs(w[e.src,k] - w[e.dst,k], MP.x[i] => val) == 0)
        end
    end
    for v in vertices(p.graph)
        @constraint(model, CP.Objective.m - subs(w[v,nₜ], t => 1) >= 0, domain = props(p.graph, v)[:cell])
    end
    @objective(model, Max, expectation([subs(w[v,1], t => 0) for v in vertices(p.graph)], μ0, p))
    return model
end


function infinite_horizon_control(CP::ControlProcess, μ0::Dict, order::Int, trange::AbstractVector, p::Partition, solver)
    MP = CP.MP
    @assert typeof(MP) == DiffusionProcess
    t = MP.time
    nₜ = length(trange)
    Δt = [trange[1], [trange[i] - trange[i-1] for i in 2:nₜ]...]
    ρ = CP.discount_factor
    model = SOSModel(solver)

    @variable(model, w[vertices(p.graph), k in 1:nₜ], Poly(monomials([MP.x..., t], 0:order)))
    @variable(model, w∞[vertices(p.graph)], Poly(monomials(MP.x, 0:order)))

    for v in vertices(p.graph), k in 1:nₜ
        X = intersect(props(p.graph, v)[:cell], @set(t >= 0 && 1-t >= 0), CP.U)
        @constraint(model, extended_inf_generator(MP, w[v,k]; scale = Δt[k]) - ρ*w[v,k] + Δt[k]*CP.Objective.l >= 0, domain = X)
    end

    for v in vertices(p.graph), k in 2:nₜ
        @constraint(model, subs(w[v,k], t => 0) - subs(w[v,k-1], t => 1) >= 0, domain=props(p.graph, v)[:cell])
    end

    for e in edges(p.graph), k in 1:nₜ
        i, val = props(p.graph, e)[:interface]
        if all(subs(MP.σ, MP.x[i] => val) .== 0)
            X = FullSpace()
            for p in props(p.graph, e.dst)[:cell].p
                if !(variables(p) == [MP.x[i]])
                    p = subs(p, MP.x[i] => val)
                    X = intersect(X, @set(p >= 0))
                end
            end
            X = intersect(X, CP.U, @set(t >= 0 && 1-t >= 0))
            @constraint(model, subs((w[e.dst,k]-w[e.src,k])*MP.f[i], MP.x[i] => val) >= 0, domain=X)
        else
            @constraint(model, subs(w[e.src,k]-w[e.dst,k], MP.x[i] => val) == 0)
        end
    end

    for v in vertices(p.graph)
        X = intersect(props(p.graph, v)[:cell], @set(t >= 0), CP.U)
        @constraint(model, extended_inf_generator(MP, w∞[v]) - ρ*w∞[v] + CP.Objective.l >= 0, domain = X)
        @constraint(model, w∞[v] - subs(w[v,nₜ], t => 1) >= 0, domain=props(p.graph, v)[:cell])
    end

    for e in edges(p.graph), k in 1:nₜ
        i, val = props(p.graph, e)[:interface]
        if all(subs(MP.σ, MP.x[i] => val) .== 0)
            X = FullSpace()
            for p in props(p.graph, e.dst)[:cell].p
                if !(variables(p) == [MP.x[i]])
                    p = subs(p, MP.x[i] => val)
                    X = intersect(X, @set(p >= 0))
                end
            end
            X = intersect(X, CP.U, @set(t >= 0))
            @constraint(model, subs((w∞[e.dst]-w∞[e.src])*MP.f[i], MP.x[i] => val) >= 0, domain=X)
        else
            @constraint(model, subs(w∞[e.src]-w∞[e.dst], MP.x[i] => val) == 0)
        end
    end
    @objective(model, Max, expectation(w, μ0, p))
    return model
end
