abstract type MarkovProcess end

mutable struct JumpProcess <: MarkovProcess
    x::Vector{<:PV} # state
    a::Vector{<:APL} # propensities
    h::Vector{Vector{<:APL}} # jumps
    X # state space enclosure
    function JumpProcess(x::Vector{<:PV}, a::Vector{<:APL}, h::Vector{Vector{T}}, X = FullSpace()) where T <: APL
        return new(x, a, h, X)
    end
end

JumpProcess(x::PV, a::APL, h::APL, X = FullSpace()) = JumpProcess([x],[a],[[h]],X)
JumpProcess(x::PV, a::Vector{<:APL}, h::Vector{<:APL}, X = FullSpace()) = JumpProcess([x],a,[[hi] for hi in h],X)

mutable struct ReactionProcess <: MarkovProcess
    ReactionSystem::ReactionSystem
    JumpProcess::JumpProcess
    species_to_index::Dict
    species_to_state::Dict
    state_to_species::Dict
    function ReactionProcess(rn::ReactionSystem)
        specs = species(rn)
        S = prodstoichmat(rn) - substoichmat(rn)
        n = length(specs)
        @polyvar(x[1:n])
        spec2idx = Dict(specs[i] => i for i in 1:n)
        spec2state = Dict(specs[i] => x[i] for i in 1:n)
        state2spec = Dict(x[i] => specs[i] for i in 1:n)
        props = reformat_reactions(reactions(rn), spec2idx, x)
        jumps = reformat_jumps(S, spec2idx, x)
        support = intersect([@set(x[i] >= 0) for i in 1:n]...)
        return new(rn, JumpProcess(x,props,jumps,support), spec2idx, spec2state, state2spec)
    end
end

mutable struct DiffusionProcess <: MarkovProcess
    x::Vector{<:PV}  # state
    f::Vector{<:APL} # drift
    σ::Matrix{<:APL} # diffusion matrix
    X # support
    function DiffusionProcess(x::Vector{<:PV}, f::Vector{<:APL}, σ::Matrix{<:APL}, X = FullSpace())
        return new(x, f, σ, X)
    end
end

DiffusionProcess(x::PV, f::APL, σ::APL, X = FullSpace()) = DiffusionProcess([x], [f], reshape([σ],1,1), X)

mutable struct JumpDiffusionProcess <: MarkovProcess
    x::Vector{<:PV} # state
    a::Vector{<:APL} # propensities
    h::Vector{Vector{<:APL}} # jumps
    f::Vector{<:APL} # drift
    σ::Matrix{<:APL} # diffusion matrix
    X # state space enclosure
    function JumpDiffusionProcess(x::Vector{<:PV}, a::Vector{<:APL}, h::Vector{Vector{<:APL}}, f::Vector{<:APL}, σ::Matrix{<:APL}, X = FullSpace())# where T <: APL
        return new(x, a, h, f, σ, X)
    end
end

JumpDiffusionProcess(x::PV, a::APL, h::APL, f::APL, σ::APL, X = Fullspace()) = JumpDiffusionProcess([x], [a], [[h]], [f], reshape(σ,1,1), X)
JumpDiffusionProcess(x::PV, a::Vector{<:APL}, h::Vector{<:APL}, f::APL, σ::APL, X = Fullspace()) = JumpDiffusionProcess([x], a, [[hi] for hi in h], [f], reshape(σ,1,1), X)

function JumpDiffusionProcess(JP::JumpProcess, DP::DiffusionProcess)
    @assert all(JP.x .== DP.x) "The jump and diffusion process must have the same state"
    return JumpDiffusionProcess(JP.x, JP.a, JP.h, DP.f, DP.σ, intersect(DP.X, JP.X))
end

mutable struct ControlProcess
    MP::MarkovProcess
    T::Real  # time horizon
    u::Vector{<:PV} # control variables
    t::PV # time variable
    U # control set
    Objective
    PathChanceConstraints
    TerminalChanceConstraints
    discount_factor
    function ControlProcess(MP, T, u, t, U, obj, PCs = [], TCs = [], dis_fac = 0)
        return new(MP, T, u, t, U, obj, PCs, TCs, dis_fac)
    end
end

ControlProcess(MP::MarkovProcess, T::Real, u::PV, t, U, obj, PCs = [], TCs = [], dis_fac = 0) =
               ControlProcess(MP, T, [u], t, U, obj, PCs, TCs, dis_fac)

mutable struct ExitProbability
    X::BasicSemialgebraicSet
end

mutable struct TerminalSetProbability
    X::BasicSemialgebraicSet
end

mutable struct LagrangeMayer
    l::APL
    m::APL
end

Lagrange(l) = LagrangeMayer(l, 0*l) # not elegant but works
Mayer(m) = LagrangeMayer(0*m, m) # not elegant but works

mutable struct ChanceConstraint
    X::BasicSemialgebraicSet
    α::Real # confidence level
end

inf_generator(MP::JumpProcess, p::Polynomial) = sum(MP.a[i]*(subs(p, MP.x => MP.h[i]) - p) for i in 1:length(MP.a))
inf_generator(MP::ReactionProcess, p::Polynomial) = inf_generator(MP.JumpProcess,p)
inf_generator(MP::DiffusionProcess, p::Polynomial) = MP.f'*∂(p,MP.x) + 1/2*sum(∂²(p,MP.x,MP.x) .* MP.σ)
inf_generator(MP::JumpDiffusionProcess, p::Polynomial) = MP.f'*∂(p,MP.x) + 1/2*sum(∂²(p,MP.x,MP.x) .* MP.σ) + sum(MP.a[i]*(subs(p, MP.x => MP.h[i]) - p) for i in 1:length(MP.a))
#old organization
#inf_generator(MP::JumpDiffusionProcess, p::Polynomial) = inf_generator(MP.JumpProcess,p) + inf_generator(MP.DiffusionProcess,p)
extended_inf_generator(MP::MarkovProcess, p::Polynomial, t::PolyVar; scale = 1) = ∂(p,t) + scale*inf_generator(MP, p)
