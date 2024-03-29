export reaction_process_setup, langevin_process_setup, stoich_bounds, transform_state, rescale_state,
       value_function, init_moments

stoichmat(rn::ReactionSystem) = prodstoichmat(rn) - substoichmat(rn)

reformat_jumps(S::Matrix, species_to_index::Dict, x::AbstractVector) = [x .+ Float64.(S[:,i]) for i in 1:size(S,2)]

∂(p,x) = differentiate.(p,x)
∂²(p,x,y) = differentiate(∂(p,x),y)
∂(X::BasicSemialgebraicSet) = [intersect(@set(X.p[i] == 0), [@set(X.p[k] >= 0) for k in 1:length(X.p) if k != i]..., X.V) for i in 1:length(X.p)]

function reformat_reactions(rxns::Vector{Reaction}, species_to_index::Dict, x, params = Dict())
    props = POLY[]
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
    if auto_scaling
        if solver isa nothing
            error("Need specification of LP solver to perform automatic scaling")
        end
        scales = stoich_bounds(rn, x0, solver)
    end
    P = ReactionProcess(rn)
    P_red = project_into_subspace(P, [x0[P.state_to_species[x]] for x in P.JumpProcess.x])
    x_scales = [scales[P_red.state_to_species[x]] for x in P_red.JumpProcess.x]
    P_final = rescale_state(P_red, x_scales)
    x0_scaled = [x0[P_red.state_to_species[x]] for x in P_red.JumpProcess.x] ./ x_scales
    return P_final, x0_scaled
end

function langevin_process_setup(rn::ReactionSystem, x0::Dict;
                                scales = Dict(s => 1.0 for s in species(rn)),
                                auto_scaling = false, solver = nothing,
                                params::Dict = Dict())
    if auto_scaling
        if solver isa nothing
            error("Need specification of LP solver to perform automatic scaling")
        end
        scales = stoich_bounds(rn, x0, solver)
    end
    P = LangevinProcess(rn)
    P_red = project_into_subspace(P, [x0[P.state_to_species[x]] for x in P.DiffusionProcess.x])
    x_scales = [scales[P_red.state_to_species[x]] for x in P_red.DiffusionProcess.x]
    P_final = rescale_state(P_red, x_scales)
    x0_scaled = [x0[P_red.state_to_species[x]] for x in P_final.DiffusionProcess.x] ./ x_scales
    return P_final, x0_scaled
end

function reaction_process_setup(rn::ReactionSystem; scales = Dict(s => 1.0 for s in species(rn)), params::Dict = Dict())
    P = ReactionProcess(rn, params)
    x_scale = Float64[scales[P.state_to_species[x]] for x in P.JumpProcess.x]
    return rescale_state(P, x_scale)
end

function langevin_process_setup(rn::ReactionSystem;
                                scales = Dict(s => 1.0 for s in species(rn)), params::Dict = Dict())
    P = LangevinProcess(rn, params)
    x_scale = [scales[P.state_to_species[x]] for x in P.DiffusionProcess.x]
    return rescale_state(P, x_scale)
end

function transform_state(P::JumpProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.x))
    Π = P.x => polynomial.(z) #[P.x[i] => z[i] for i in 1:length(P.x)]
    a = map(p -> subs(p,Π), P.a)
    h = [map(p -> subs(p,Π), h[iv]) for h in P.h]
    X = subs_X(P.X, Π) #P.X isa FullSpace ? P.X : intersect([@set(p(Π) >= 0) for p in P.X.p]...)
    return JumpProcess(x, a, h, X, P.iv, P.controls, P.poly_vars)
end

function transform_state(P::ReactionProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.JumpProcess.x))
    Π = P.JumpProcess.x => polynomial.(z) #[P.JumpProcess.x[i] => z[i] for i in 1:length(P.JumpProcess.x)]
    species_to_state = Dict(spec => P.species_to_state[spec](Π) for spec in keys(P.species_to_state))
    #state_to_species = Dict(P.species_to_state[spec] => spec for spec in keys(P.species_to_state)) # this is wrong, need to update with the new independent variabls x
    state_to_species = Dict(xi => P.state_to_species[xi] for xi in x)
    return ReactionProcess(P.ReactionSystem, transform_state(P.JumpProcess, x, z; iv = iv), P.species_to_index, species_to_state, state_to_species)
end

function transform_state(P::DriftProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.x))
    Π = P.x => polynomial.(z)
    f = map(p -> subs(p,Π), P.f[iv])
    X = subs_X(P.X, Π) #P.X isa FullSpace ? P.X : intersect([@set(p(Π) >= 0) for p in P.X.p]...)
    return DriftProcess(x, f, X, P.iv, P.controls, P.poly_vars)
end

function transform_state(P::DiffusionProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.x))
    Π = P.x => polynomial.(z)
    f = map(p -> subs(p,Π), P.f[iv])
    σ = map(p -> subs(p,Π), P.σ[iv, iv])
    X = subs_X(P.X, Π) #P.X isa FullSpace ? P.X : intersect([@set(p(Π) >= 0) for p in P.X.p]...)
    return DiffusionProcess(x, f, σ, X, P.iv, P.controls, P.poly_vars)
end

function transform_state(P::LangevinProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.DiffusionProcess.x))
    Π = P.DiffusionProcess.x => polynomial.(z)
    species_to_state = Dict(spec => P.species_to_state[spec](Π) for spec in keys(P.species_to_state))
    state_to_species = Dict(P.species_to_state[spec] => spec for spec in keys(P.species_to_state))
    diff_process = transform_state(P.DiffusionProcess, x, z; iv = iv)
    return LangevinProcess(P.ReactionSystem, diff_process, P.species_to_index, species_to_state, state_to_species)
end

function transform_state(P::JumpDiffusionProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.x))
    Π = P.x => polynomial.(z) #[P.x[i] => z[i] for i in 1:length(P.x)]
    a = map(p -> subs(p,Π), P.a)
    h = [map(p -> subs(p,Π), h[iv]) for h in P.h]
    f = map(p -> subs(p,Π), P.f[iv])
    σ = map(p -> subs(p,Π), P.σ[iv, iv])
    X = subs_X(P.X, Π) #P.X isa FullSpace ? P.X : intersect([@set(p(Π) >= 0) for p in P.X.p]...)
    return JumpDiffusionProcess(x, a, h, f, σ, X, P.iv, P.controls, P.poly_vars)
end

function rescale_state(P::JumpProcess, x0::Vector{<:Real})
    rescaled_P = transform_state(P, P.x, P.x .* x0)
    for i in 1:length(new_P.h)
        rescaled_P.h[i] ./= x0
    end
    return rescaled_P
end

function rescale_state(P::ReactionProcess, x0::Vector{<:Real})
    rescaled_P = transform_state(P, P.JumpProcess.x, polynomial.(P.JumpProcess.x .* x0))
    for i in 1:length(P.JumpProcess.h)
        rescaled_P.JumpProcess.h[i] ./= x0
    end
    return rescaled_P
end

function rescale_state(P::DriftProcess, x0::Vector{<:Real})
    rescaled_P = transform_state(P, P.x, P.x .* x0)
    rescaled_P.f ./= x0
    return rescaled_P
end

function rescale_state(P::DiffusionProcess, x0::Vector{<:Real})
    rescaled_P = transform_state(P, P.x, P.x .* x0)
    rescaled_P.f ./= x0
    rescaled_P.σ ./= x0*x0'
    return rescaled_P
end

function rescale_state(P::LangevinProcess, x0::Vector{<:Real})
    rescaled_P = transform_state(P, P.DiffusionProcess.x, polynomial.(P.DiffusionProcess.x .* x0))
    rescaled_P.DiffusionProcess.f ./= x0
    rescaled_P.DiffusionProcess.σ ./= x0*x0'
    return rescaled_P
end

function rescale_state(P::JumpDiffusionProcess, x0::Vector{<:Real})
    rescaled_P = transform_state(P, P.x, P.x .* x0)
    for i in 1:length(P.h)
        rescaled_P.h[i] ./= x0
    end
    rescaled_P.f ./= x0
    rescaled_P.σ ./= x0*x0'
    return rescaled_P
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

function project_into_subspace(P::MarkovProcess, B, f::Vector{<:Real})
    iv, dv, Q, R = partition_variables(B)
    if !isempty(dv)
        z = Vector{POLY}(undef, size(B,2))
        x = states(P)[iv] #@polyvar(x[1:length(iv)])
        z[iv] .= polynomial.(x)
        z[dv] .= polynomial.(round.(R[:,dv]\(Q'*f), digits = 12) - round.((R[:,dv]\R[:,iv]), digits = 12)*x)
        projected_P = transform_state(P, x, z; iv=iv)
    else
        projected_P = P
    end
    return projected_P
end

function project_into_subspace(P::Union{ReactionProcess,LangevinProcess}, x0::Vector{<:Real})
    B = transpose(nullspace(transpose(stoichmat(P.ReactionSystem))))
    if !isempty(B)
        f = B*x0
        projected_P = project_into_subspace(P, B, f)
    else
        projected_P = P
    end
    return projected_P
end

function init_moments(x, x0, d)
    return Dict(mono => mono(x=>x0) for mono in  monomials(x, 0:d))
end

init_moments(x, x0, d, P::Partition) =  Dict(P.get_vertex(x0) => init_moments(x, x0, d))

function expectation(w::APL, μ)
    ex = AffExpr(0.0)
    for i in 1:length(w.a)
        if w.a[i] != 0 && μ[w.x[i]] != 0
            ex += w.a[i] * μ[w.x[i]]
        end
    end
    return ex
end

expectation(w::VariableRef, μ) = w*μ[1]
expectation(w::AbstractVector, μ) = sum(expectation(w[v], μ[v]) for v in keys(μ))
expectation(w::Dict, μ) = sum(expectation(w[v], μ[v]) for v in keys(μ))

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
    eqs = POLY[]
    for eq in equalities(X)
        push!(eqs, eq(submap))
    end
    ineqs = POLY[]
    for eq in inequalities(X)
        push!(ineqs, eq(submap))
    end
    return BasicSemialgebraicSet(algebraic_set(eqs), ineqs)
end

function subs_X(X::FullSpace, submap)
    return X
end

# applies only for constant jumps
function reverse_jumps(jumps)
    rev_jumps = deepcopy(jumps)
    for j in rev_jumps
        for e in j 
            @assert (maxdegree(e) <= 1 && e.a[end] == 1) "jump reversal currently only supported for constant jumps"
            e.a[1] = length(e.a) > 1 ? -1*e.a[1] : e.a[1]
        end
    end
    return rev_jumps
end

inequalities(::FullSpace) = []
polynomial(a::Real) = polynomial(DynPoly.Term(a, Monomial(1))) #piracy ... do we need this?

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

struct StateDict{dType,kType}
    dict::dType
    keys::kType
    tol::Float64
end

function StateDict(entries::AbstractVector{<:Pair}; tol=1e-8)
    dict = Dict(entries)
    return StateDict(dict, collect(keys(dict)), convert(Float64, tol))
end

function StateDict((entries::Pair)...; tol = 1e-8)
    dict = Dict(entries)
    return StateDict(dict, collect(keys(dict)), convert(Float64, tol))
end

import Base.getindex, Base.setindex!, Base.values, Base.keys

function getindex(sd::StateDict, idx)
    sd_idx = findfirst(key -> isapprox(key,idx,atol=sd.tol), sd.keys)
    return sd.dict[sd.keys[sd_idx]] 
end

function setindex!(sd::StateDict, val, idx)
    sd.dict[idx] = val
    push!(sd.keys, idx)
end

values(sd::StateDict) = values(sd.dict)
keys(sd::StateDict) = keys(sd.dict)

function closeto(x, sd::StateDict)
    sd_idx = findfirst(key -> isapprox(key,x,atol=sd.tol), sd.keys)
    return !isnothing(sd_idx)
end

states(P::MarkovProcess) = P.x
states(P::ReactionProcess) = P.JumpProcess.x
states(P::LangevinProcess) = P.DiffusionProcess.x
