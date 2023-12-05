export grid_graph, grid_partition, discrete_grid_graph, discrete_grid_partition, props

function grid_graph(x, lb, ub, n; inf_top = zeros(Int64, length(ub)), inf_floor = zeros(Int64, length(lb)), X_base = FullSpace())
    @assert all(lb .<= ub) "lower boundary must be smaller than upper boundary"
    @assert all(n .>= 1) "minimum one partition element per dimension required"
    x_ranges = []
    for k in 1:length(n)
        n_eff = n[k] - inf_top[k] - inf_floor[k] # number of partition elements between lb[k] and ub[k]
        if n_eff < 0
            @warn "number of partition elements in dimension $k is inconsistent with (semi-)infinite domain; dimension ignored for partitioning!"
            push!(x_ranges, [inf_floor[k] == 0 ? lb[k] : -Inf, (inf_top[k] == 0 ? ub[k] : Inf)])
        elseif n_eff == 0 
            if inf_floor[k] + inf_top[k] == 1
                push!(x_ranges, [inf_floor[k] == 0 ? lb[k] : -Inf, (inf_top[k] == 0 ? ub[k] : Inf)])
            else 
                if lb[k] != ub[k] 
                    @warn "If 2 partition elements are specified on infinite domain, \n resolved range [a,b] must be a singleton (a=b) for consistency. The partition then reads [-Inf, a, Inf]. \n 
                           Here the resolved range does not satisfy this requirement. Proceed with partition [-Inf, (a+b)/2, Inf]."
                end
                mid_point = (lb[k] + ub[k])/2
                push!(x_ranges, [-Inf, mid_point, Inf])
            end
        else
            boundaries = Float64[]
            if inf_floor[k] == 1
                push!(boundaries, -Inf)
            end
            boundaries = vcat(boundaries, collect(range(lb[k], stop=ub[k], length=n_eff+1)))
            if inf_top[k] == 1
                push!(boundaries, Inf)
            end
            push!(x_ranges, boundaries)
        end
    end
    return grid_graph(x, x_ranges; X_base = X_base)
end

function grid_graph(x, x_ranges; X_base = FullSpace())
    n = length.(x_ranges) .- 1
    nv = prod(n)

    mg = MetaDiGraph(nv)
    for i in 1:nv
        idx = invert_index(i, n)
        subset = X_base
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
                set_prop!(mg, Edge(i,j), :interface, [intersect(@set(x[k] == x_ranges[k][idx[k]]), X_base)])
                idx[k] -= 1
            end
        end
    end
    get_vertex = function (x)
                    idx = zeros(Int64, length(n))
                    in_set = true
                    for k in 1:length(n)
                        j = findfirst(m -> m >= x[k], x_ranges[k]) - 1
                        if isnothing(j)
                            @warn "state outside state space"
                            in_set = false
                            break
                        end
                        idx[k] = max(j, 1)
                    end
                    vertex = in_set ? linearize_index(idx, n) : -1

                    return vertex
                 end
    return mg, get_vertex
end


function grid_partition(x, lb, ub, n; inf_top = zeros(Int64, length(ub)), inf_floor = zeros(Int64, length(lb)), X_base = FullSpace())
    return Partition(grid_graph(x, lb, ub, n; inf_top = inf_top, inf_floor = inf_floor, X_base = X_base)...)
end

function grid_partition(x, x_ranges, X_base = FullSpace())
    return Partition(grid_graph(x, x_ranges, X_base = X_base)...)
end

function check_membership(X::AbstractSemialgebraicSet, x_var, x)
    is_member = check_membership(X.V, x_var, x)
    if is_member
        for ineq in inequalities(X)
            if ineq(x_var => x) < 0
                is_member = false
                break
            end
        end
    end
    return is_member
end

function check_membership(X::AbstractAlgebraicSet, x_var, x)
    is_member = true
    for eq in equalities(X)
        if eq(x_var => x) != 0
            is_member = false
            break
        end
    end
    return is_member
end

check_membership(X::FullSpace, x) = true

# not elegant but at least correct!
function complement(X::BasicSemialgebraicSet, H = FullSpace())
    Ys = BasicSemialgebraicSet[]
    @assert isempty(equalities(X)) "Equality constraints are not allowed when computing the complement"
    ineqs = inequalities(X)
    m = length(ineqs)
    for i in 1:m
        push!(Ys, intersect(H, BasicSemialgebraicSet(algebraic_set(POLY[]),
                                                     vcat(ineqs[1:end-1], -ineqs[end]))))
        pop!(ineqs)
    end
    return Ys
end

function discrete_grid_graph(JP::JumpProcess, x_ranges, Xc)
    x_var = JP.x
    partition_graph = MetaDiGraph()
    discrete_states = collect(product(x_ranges...))
    state_to_vertex = StateDict(;tol=1e-3)
    n_discrete = length(discrete_states)
    for i in eachindex(discrete_states)
        add_vertex!(partition_graph, :cell, @singleton(discrete_states[i]))
        state_to_vertex[[discrete_states[i]...]] = i 
    end
    add_vertex!(partition_graph, :cell, Xc)

    rev_jumps = reverse_jumps(JP.h)
    for i in eachindex(discrete_states)
        neighbors = neighborhood(discrete_states[i], JP, rev_jumps, state_to_vertex)
        if !isempty(neighbors)
            add_edge!(partition_graph, n_discrete + 1, state_to_vertex[[discrete_states[i]...]])
        end
    end

    get_vertex = function (x)
        if closeto(x, state_to_vertex)
            vertex = state_to_vertex[x]
        elseif check_membership(JP.X, x_var, x)
            vertex = n_discrete + 1 
        else
            vertex = -1 
        end
        return vertex
    end

    return partition_graph, get_vertex
end

function discrete_grid_partition(JP::JumpProcess, x_ranges, Xc)
    return Partition(discrete_grid_graph(JP, x_ranges, Xc)...)
end

function neighborhood(state, JP::JumpProcess, state_to_vertex)
    x_var = JP.x
    rev_jumps = reverse_jumps(JP.h)
    neighbors = []
    for i in eachindex(rev_jumps)
        origin = [jump_component(x_var => state) for jump_component in rev_jumps[i]]
        if subs(JP.a[i], x_var => origin) != 0 && !closeto(origin, state_to_vertex) && check_membership(JP.X, x_var, origin)
            push!(neighbors, origin)
        end
    end
    return neighbors
end

function neighborhood(state, JP::JumpProcess, rev_jumps, state_to_vertex)
    x_var = JP.x
    neighbors = []
    for i in eachindex(rev_jumps)
        origin = round.([jump_component(x_var => state) for jump_component in rev_jumps[i]], digits = 8)
        if subs(JP.a[i], x_var => origin) != 0 && !closeto(origin, state_to_vertex) && check_membership(JP.X, x_var, origin)
            push!(neighbors, origin)
        end
    end
    return neighbors
end