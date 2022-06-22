export reaction_process_setup, langevin_process_setup, stoich_bounds, transform_state!, rescale_state!,
       value_function

stoichmat(rn::ReactionSystem) = prodstoichmat(rn) - substoichmat(rn)

reformat_jumps(S::Matrix, species_to_index::Dict, x::AbstractVector) = [x .+ S[:,i] for i in 1:size(S,2)]

∂(p,x) = differentiate.(p,x)
∂²(p,x,y) = differentiate(∂(p,x),y)
∂(X::BasicSemialgebraicSet) = [intersect(@set(X.p[i] == 0), [@set(X.p[k] >= 0) for k in 1:length(X.p) if k != i]..., X.V) for i in 1:length(X.p)]

function reformat_reactions(rxns::Vector{Reaction}, species_to_index::Dict, x::Vector{<:PV}, params = Dict())
    props = Polynomial{true,Float64}[]
    for r in rxns
        @unpack rate, substrates, substoich, only_use_rate = r
        rate = rate isa Sym ? params[rate] : rate
        a = rate*polynomial(MonomialVector(x,0))
        if !only_use_rate
            for (s, ν) in enumerate(substoich)
                idx = species_to_index[substrates[s]]
                a *= prod(x[idx] - i for i in 0:ν-1)/factorial(ν) # consistent with rxns
            end
        end
        push!(props, a)
    end
    return props
end

function stoich_bounds(rn::ReactionSystem, x0::Dict, solver)
    S = stoichmat(rn)
    B = nullspace(S)
    specs = species(rn)
    x_init = [x0[spec] for spec in specs]
    m = Model(solver)
    @variable(m, x[1:length(specs)] >= 0)
    @constraint(m, B'*x .== B'*x_init)
    scales = Dict()
    for i in 1:length(specs)
        @objective(m, Max, x[i])
        optimize!(m)
        if termination_status(m) == MOI.DUAL_INFEASIBLE
            @warn  "state space unbounded, autoscaling failed!"
            scales[specs[i]] = 1.0
        else
            scales[specs[i]] = min(1.0, objective_value(m))
        end
    end
    return scales
end

"""
    reaction_process_setup(rn::ReactionSystem, x0::Dict;
                    scales = Dict(s => 1.0 for s in species(rn)),
                    auto_scaling = false, solver = nothing)

    transforms a reaction network as defined by Catalyst.jl's ReactionSystem type
    into an equivalent ReactionProcess accounting for reaction invariants and
    scales.
"""
function reaction_process_setup(rn::ReactionSystem, x0::Dict;
                                scales = Dict(s => 1.0 for s in species(rn)),
                                auto_scaling = false, solver = nothing,
                                params::Dict = Dict())
    P = ReactionProcess(rn, params)
    project_into_subspace!(P, Float64[x0[s] for s in species(rn)])
    if auto_scaling
        if solver isa nothing
            error("Need specification of LP solver to perform automatic scaling")
        end
        scales = stoich_bounds(rn, x0, solver)
    end
    x_scale = [scales[P.state_to_species[x]] for x in P.JumpProcess.x]
    x0 = [x0[P.state_to_species[x]]/scales[P.state_to_species[x]] for x in P.JumpProcess.x]
    rescale_state!(P, x_scale)
    return P, x0
end

function langevin_process_setup(rn::ReactionSystem, x0::Dict;
                                scales = Dict(s => 1.0 for s in species(rn)),
                                auto_scaling = false, solver = nothing,
                                params::Dict = Dict())
    P = LangevinProcess(rn, params)
    project_into_subspace!(P, Float64[x0[s] for s in species(rn)])
    if auto_scaling
        if solver isa nothing
            error("Need specification of LP solver to perform automatic scaling")
        end
        scales = stoich_bounds(rn, x0, solver)
    end
    x_scale = [scales[P.state_to_species[x]] for x in P.DiffusionProcess.x]
    x0 = [x0[P.state_to_species[x]]/scales[P.state_to_species[x]] for x in P.DiffusionProcess.x]
    rescale_state!(P, x_scale)
    return P, x0
end

function reaction_process_setup(rn::ReactionSystem; scales = Dict(s => 1.0 for s in species(rn)), params::Dict = Dict())
    P = ReactionProcess(rn, params)
    x_scale = Float64[scales[P.state_to_species[x]] for x in P.JumpProcess.x]
    rescale_state!(P, x_scale)
    return P
end

function langevin_process_setup(rn::ReactionSystem;
                                scales = Dict(s => 1.0 for s in species(rn)), params::Dict = Dict())
    P = LangevinProcess(rn, params)
    x_scale = [scales[P.state_to_species[x]] for x in P.DiffusionProcess.x]
    rescale_state!(P, x_scale)
    return P
end

function transform_state!(P::JumpProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.x))
    Π = P.x => polynomial.(z) #[P.x[i] => z[i] for i in 1:length(P.x)]
    P.a = map(p -> p(Π), P.a)
    P.h = [map(p -> p(Π), h[iv]) for h in P.h]
    P.X = subs_X(P.X, Π) #P.X isa FullSpace ? P.X : intersect([@set(p(Π) >= 0) for p in P.X.p]...)
    P.x = x
end

function transform_state!(P::ReactionProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.JumpProcess.x))
    Π = P.JumpProcess.x => polynomial.(z) #[P.JumpProcess.x[i] => z[i] for i in 1:length(P.JumpProcess.x)]
    P.species_to_state = Dict(spec => P.species_to_state[spec](Π) for spec in keys(P.species_to_state))
    P.state_to_species = Dict(P.species_to_state[spec] => spec for spec in keys(P.species_to_state))
    transform_state!(P.JumpProcess, x, z; iv = iv)
end

function transform_state!(P::DiffusionProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.x))
    Π = P.x => polynomial.(z)
    P.f = map(p -> p(Π), P.f[iv])
    P.σ = map(p -> p(Π), P.σ[iv, iv])
    P.X = subs_X(P.X, Π) #P.X isa FullSpace ? P.X : intersect([@set(p(Π) >= 0) for p in P.X.p]...)
    P.x = x
end

function transform_state!(P::LangevinProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.DiffusionProcess.x))
    Π = P.DiffusionProcess.x => polynomial.(z)
    P.species_to_state = Dict(spec => P.species_to_state[spec](Π) for spec in keys(P.species_to_state))
    P.state_to_species = Dict(P.species_to_state[spec] => spec for spec in keys(P.species_to_state))
    transform_state!(P.DiffusionProcess, x, z; iv = iv)
end

function transform_state!(P::JumpDiffusionProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.x))
    Π = P.x => polynomial.(z) #[P.x[i] => z[i] for i in 1:length(P.x)]
    P.a = map(p -> p(Π), P.a)
    P.h = [map(p -> p(Π), h[iv]) for h in P.h]
    P.f = map(p -> p(Π), P.f[iv])
    P.σ = map(p -> p(Π), P.σ[iv, iv])
    P.X = subs_X(P.X, Π) #P.X isa FullSpace ? P.X : intersect([@set(p(Π) >= 0) for p in P.X.p]...)
    P.x = x
end

function rescale_state!(P::JumpProcess, x0::Vector{<:Real})
    transform_state!(P, P.x, P.x .* x0)
    for i in 1:length(P.h)
        P.h[i] ./= x0
    end
end

function rescale_state!(P::ReactionProcess, x0::Vector{<:Real})
    transform_state!(P, P.JumpProcess.x, polynomial.(P.JumpProcess.x .* x0))
    for i in 1:length(P.JumpProcess.h)
        P.JumpProcess.h[i] ./= x0
    end
end

function rescale_state!(P::DiffusionProcess, x0::Vector{<:Real})
    transform_state!(P, P.x, P.x .* x0)
    P.f ./= x0
    P.σ ./= x0*x0'
end

function rescale_state!(P::LangevinProcess, x0::Vector{<:Real})
    transform_state!(P, P.DiffusionProcess.x, polynomial.(P.DiffusionProcess.x .* x0))
    P.DiffusionProcess.f ./= x0
    P.DiffusionProcess.σ ./= x0*x0'
end

function rescale_state!(P::JumpDiffusionProcess, x0::Vector{<:Real})
    transform_state!(P, P.x, P.x .* x0)
    for i in 1:length(P.h)
        P.h[i] ./= x0
    end
    P.f ./= x0
    P.σ ./= x0*x0'
end

function partition_variables(B)
    m, n = size(B)
    @assert(m <= n)
    Q, R = qr(B)
    ds = []
    i, k = 1, 1
    while k <= n && i <= m
        if abs(R[i,k]) >= 1e-12
            push!(ds, k)
            i += 1
        end
        k += 1
    end
    is = setdiff(1:n, ds)
    return is, ds, Q, R
end

function project_into_subspace!(P::MarkovProcess, B, f::Vector{<:Real})
    iv, dv, Q, R = partition_variables(B)
    if !isempty(dv)
        z = Vector{Polynomial{true, Float64}}(undef, size(B,2))
        @polyvar(x[1:length(iv)])
        z[iv] .= polynomial.(x)
        z[dv] .= polynomial.(round.(R[:,dv]\(Q'*f), digits = 12) - round.((R[:,dv]\R[:,iv]), digits = 12)*x)
        transform_state!(P, x, z; iv=iv)
    end
end

function project_into_subspace!(P::Union{ReactionProcess,LangevinProcess}, x0::Vector{<:Real})
    B = transpose(nullspace(transpose(stoichmat(P.ReactionSystem))))
    if !isempty(B)
        f = B*x0
        project_into_subspace!(P, B, f)
    end
end

function init_moments(x, x0, d)
    n = length(x0)
    moms = Dict()
    for i in 1:(d+1)^n
        midx = invert_index(i, (d+1)*ones(Int64,n)) .- 1
        moms[prod(x.^midx)] = prod(x0.^midx)
    end
    return moms
end

init_moments(x, x0, d, P::Partition) =  Dict(P.get_vertex(x0) => init_moments(x, x0, d))

linearize_index(idx,rs) = idx[1] + (length(idx) > 1 ? sum((idx[i] - 1) * prod(rs[1:i-1]) for i in 2:length(idx)) : 0)

function invert_index(idx, rs)
    inv_idx = similar(rs)
    for i in length(rs):-1:2
        fac = prod(rs[1:i-1])
        n = div(idx - 1, fac)
        inv_idx[i] = n + 1
        idx -= n*fac
    end
    inv_idx[1] = idx
    return inv_idx
end

function expectation(w::Polynomial, μ::Dict)
    ex = AffExpr(0.0)
    for i in 1:length(w.a)
        if w.a[i] != 0 && μ[w.x[i]] != 0
            ex += w.a[i] * μ[w.x[i]]
        end
    end
    return ex
end

expectation(w::VariableRef, μ::Dict) = w*μ[1]
expectation(w::AbstractVector, μ::Dict) = sum(expectation(w[v], μ[v]) for v in keys(μ))
expectation(w::Dict, μ::Dict) = sum(expectation(w[v], μ[v]) for v in keys(μ))

function value_function(CP::ControlProcess, trange, bound::Bound)
    if trange[1] == 0
        trange = trange[2:end]
    end
    V_poly(t,x) = t <= trange[end] ? bound.w[get_piece(t, trange), bound.partition.get_vertex(x)] : bound.w[bound.partition.get_vertex(x)]
    V_val(t,x) = V_poly(t,x)(CP.MP.iv => t, CP.MP.x => x)
    return V_poly, V_val
end

function get_piece(t, trange)
    if t > trange[end] || t < 0
        @warn string("t = ", t, " outside the domain of V")
        i = t < 0 ? 1 : length(trange)
    elseif t == trange[end]
        i = length(trange)
    else
        i = findfirst(ti -> ti > t, trange)
    end
    return i
end

function trivial_partition(X::AbstractSemialgebraicSet)
    G = MetaDiGraph()
    add_vertex!(G, :cell, X)
    return Partition(G, x -> 1)
end

function split_state_space(MP::MarkovProcess, X::BasicSemialgebraicSet)
    Xc = complement(X, MP.X)
    G = MetaDiGraph()
    add_vertex!(G, :cell, X)
    add_vertex!(G, :cell, Xc)
    add_edge!(G, 1, 2, :interface, ∂(X))
    return Partition(G, x -> check_membership(x, X) ? 1 : 2)
end

function dual_poly(w, t, trange)
    if trange[1] == 0
        trange = trange[2:end]
    end
    return Dict(key => subs(value(w[key]), t => key[1] > 1 ? (t - trange[key[1]-1])/(trange[key[1]] - trange[key[1]-1]) : t/trange[key[1]]) for key in keys(w))
end

function dual_poly(w)
    return Dict(key => value(w[key]) for key in keys(w))
end

function subs_X(X::BasicSemialgebraicSet, submap)
    eqs = Polynomial{true, Float64}[]
    for eq in equalities(X)
        push!(eqs, eq(submap))
    end
    ineqs = Polynomial{true, Float64}[]
    for eq in inequalities(X)
        push!(ineqs, eq(submap))
    end
    return BasicSemialgebraicSet(algebraicset(eqs), ineqs)
end

inequalities(::FullSpace) = []
polynomial(a::Real) = polynomial(DPTerm{true}(a))