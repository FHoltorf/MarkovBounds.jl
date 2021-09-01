export grid_graph, Partition, props, Singleton, @singleton

mutable struct Singleton <: AbstractSemialgebraicSet
    x::AbstractVector
end

Singleton(x::Real) = Singleton([x])
Singleton(x::Tuple) = Singleton([x...])

macro singleton(x)
    return :(Singleton($(esc(x))))
end

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
                set_prop!(mg, Edge(i,j), :interface, [@set(x[k] == x_ranges[k][idx[k]])])
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

function check_membership(X::AbstractSemialgebraicSet, x_var, x)
    is_member = check_membership(X.V)
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

# not elegant but at least in principle correct!
function complement(X::BasicSemialgebraicSet, H = FullSpace())
    Ys = BasicSemialgebraicSet[]
    @assert isempty(equalities(X)) "Equality constraints are not allowed when computing the complement"
    ineqs = inequalities(X)
    m = length(ineqs)
    for i in 1:m
        push!(Ys, intersect(H, BasicSemialgebraicSet(algebraicset(Polynomial{true, Float64}[]),
                                                    vcat(ineqs[1:end-1], -ineqs[end]))))
        pop!(ineqs)
    end
    return Ys
end

function subs_X(X::BasicSemialgebraicSet, submap)
    eqs = Polynomial{true, Float64}[]
    for eq in equalities(X)
        push!(eqs, subs(eq, submap))
    end
    ineqs = Polynomial{true, Float64}[]
    for eq in inequalities(X)
        push!(ineqs, subs(eq, submap))
    end
    return BasicSemialgebraicSet(algebraicset(eqs), ineqs)
end

# in general this is conservative but holds
