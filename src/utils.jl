stoichmat(rn::ReactionSystem) = prodstoichmat(rn) - substoichmat(rn)

reformat_jumps(S::Matrix, species_to_index::Dict, x::AbstractVector) = [x .+ S[i,:] for i in 1:size(S,1)]

∂(p,x) = differentiate.(p,x)
∂²(p,x,y) = differentiate(∂(p,x),y)
∂(X) = [intersect(@set(X.p[i] == 0), [@set(X.p[k] >= 0) for k in 1:length(X.p) if k != i]..., X.V) for i in 1:length(X.p)]

function reformat_reactions(rxns::Vector{Reaction}, species_to_index::Dict, x::Vector{<:PV})
    props = APL[]
    for r in rxns
        @unpack rate, substrates, substoich, only_use_rate = r
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

function setup_reaction_process(rn::ReactionSystem, x0::Dict; scales = Dict(s => 1.0 for s in species(rn)), auto_scaling = false, solver = nothing)
    P = ReactionProcess(rn)
    project_into_subspace!(P, [x0[s] for s in species(rn)])
    if auto_scaling
        scales = stoich_bounds(rn, x0, solver)
    end
    x_scale = [scales[P.state_to_species[x]] for x in P.JumpProcess.x]
    x0 = [x0[P.state_to_species[x]]/scales[P.state_to_species[x]] for x in P.JumpProcess.x]
    rescale_state!(P, x_scale)
    return P, x0
end

function setup_reaction_process(rn::ReactionSystem; scales = Dict(s => 1.0 for s in species(rn)))
    P = ReactionProcess(rn)
    x_scale = [scales[P.state_to_species[x]] for x in P.JumpProcess.x]
    rescale_state!(P, x_scale)
    return P
end

function transform_state!(P::JumpProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.x))
    Π = [P.x[i] => z[i] for i in 1:length(P.x)]
    P.a = subs.(P.a, Π...)
    P.h = [subs.(h[iv], Π...) for h in P.h]
    P.X = intersect([@set(subs.(p, Π...) >= 0) for p in P.X.p]...)
    P.x = x
end

function transform_state!(P::ReactionProcess, x::AbstractVector, z::AbstractVector; iv::AbstractVector = 1:length(P.JumpProcess.x))
    Π = [P.JumpProcess.x[i] => z[i] for i in 1:length(P.JumpProcess.x)]
    P.species_to_state = Dict(spec => subs(P.species_to_state[spec], Π...) for spec in keys(P.species_to_state))
    P.state_to_species = Dict(P.species_to_state[spec] => spec for spec in keys(P.species_to_state))
    transform_state!(P.JumpProcess, x, z; iv = iv)
end

function rescale_state!(P::DiffusionProcess, x0::Vector{<:Number})
    transform_state!(P, P.x, P.x .* x0)
    P.f ./= x0
    P.σ ./= x0*x0'
end

function rescale_state!(P::JumpProcess, x0::Vector{<:Number})
    transform_state!(P, P.x, P.x .* x0)
    for i in 1:length(P.h)
        P.h[i] ./= x0
    end
end

function rescale_state!(P::JumpDiffusionProcess, x0::Vector{<:Number})
    rescale_state!(P.JumpProcess.x, x0)
    rescale_state!(P.DiffusionProcess.x, x0)
end

function rescale_state!(P::ReactionProcess, x0::Vector{<:Number})
    transform_state!(P, P.JumpProcess.x, P.JumpProcess.x .* x0)
    for i in 1:length(P.JumpProcess.h)
        P.JumpProcess.h[i] ./= x0
    end
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

function project_into_subspace!(P::MarkovProcess, B, f::Vector{<:Number})
    iv, dv, Q, R = partition_variables(B)
    if !isempty(dv)
        z = Vector{Polynomial}(undef, size(B,2))
        @polyvar(x[1:length(iv)])
        z[iv] .= polynomial.(x)
        z[dv] .= polynomial.(round.(R[:,dv]\(Q'*f), digits = 12) - round.((R[:,dv]\R[:,iv]), digits = 12)*x)
        transform_state!(P, x, z; iv=iv)
    end
end

function project_into_subspace!(P::ReactionProcess, x0::Vector{<:Number})
    B = transpose(nullspace(stoichmat(P.ReactionSystem)))
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

function init_moments(x, x0, d, p::Partition)
    v0 = p.get_vertex(x0)
    n = length(x0)
    moms = Dict()
    for i in 1:(d+1)^n
        midx = invert_index(i, (d+1)*ones(Int64,n)) .- 1
        for v in vertices(p.graph)
            moms[v,prod(x.^midx)] = (v == v0 ? prod(x0.^midx) : 0.0)
        end
    end
    return moms
end

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
        if !(w.a[i] == 0)
            ex += w.a[i] * μ[w.x[i]]
        end
    end
    return ex
end

function expectation(w, μ::Dict, p::Partition)
    ex = AffExpr(0.0)
    for v in vertices(p.graph)
        for i in 1:length(w[v,1].a)
            if μ[v, w[v].x[i]] != 0
                ex += w[v].a[i]*μ[v, w[v].x[i]]
            end
        end
    end
    return ex
end
