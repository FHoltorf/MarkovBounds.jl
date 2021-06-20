abstract type MarkovProcess end

mutable struct JumpProcess <: MarkovProcess
    x::Vector{<:PV} # state
    a::Vector{<:APL} # propensities
    h::Vector{Vector{<:APL}} # jumps
    X # state space enclosure
    time::PV
    controls::Vector{<:PV}
    poly_vars::Dict # needed if process defined in terms of symbolics.jl variables
    function JumpProcess(x::Vector{<:PV}, a::Vector{<:APL}, h::Vector{<:Vector{<:APL}}, X = FullSpace();
                         time = @polyvar(t)[1], controls = PV{true}[], poly_vars = Dict())
         return new(x, a, h, X, time, controls, poly_vars)
    end
end

JumpProcess(x::PV, a::APL, h::APL, X = FullSpace(); time = @polyvar(t)[1], controls = PV{true}[]) =
            JumpProcess([x], [a], [[h]], X; time=time, controls = (controls isa Vector ? controls : [controls]))
JumpProcess(x::PV, a::Vector{<:APL}, h::Vector{<:APL}, X = FullSpace(); time = @polyvar(t)[1], controls = PV{true}[]) =
            JumpProcess([x], a, [[hi] for hi in h], X; time = time, controls = (controls isa Vector ? controls : [controls]))

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
    time::PV
    controls::Vector{<:PV}
    poly_vars::Dict # needed if process defined in terms of symbolics.jl variables
    function DiffusionProcess(x::Vector{<:PV}, f::Vector{<:APL}, σ::Matrix{<:APL}, X = FullSpace();
                              time = @polyvar(t)[1], controls = PV{true}[], poly_vars = Dict())
        return new(x, f, σ, X, time, controls, poly_vars)
    end
end

DiffusionProcess(x::PV, f::APL, σ::APL, X = FullSpace(); time::PV = @polyvar(t)[1], controls = PV{true}[]) =
                 DiffusionProcess([x], [f], reshape([σ],1,1), X; time=time, controls = (controls isa Vector ? controls : [controls]))

mutable struct JumpDiffusionProcess <: MarkovProcess
    x::Vector{<:PV} # state
    a::Vector{<:APL} # propensities
    h::Vector{Vector{<:APL}} # jumps
    f::Vector{<:APL} # drift
    σ::Matrix{<:APL} # diffusion matrix
    X # state space enclosure
    time::PV
    controls::Vector{<:PV}
    poly_vars::Dict # needed if process defined in terms of symbolics.jl variables
    function JumpDiffusionProcess(x::Vector{<:PV}, a::Vector{<:APL}, h::Vector{<:Vector{<:APL}}, f::Vector{<:APL}, σ::Matrix{<:APL}, X = FullSpace();
                                  time = @polyvar(t)[1], controls = PV{true}[], poly_vars = Dict())
        return new(x, a, h, f, σ, X, time, controls, poly_vars)
    end
end

JumpDiffusionProcess(x::PV, a::APL, h::APL, f::APL, σ::APL, X = Fullspace(); time = @polyvar(t)[1], controls = PV{true}[]) =
                     JumpDiffusionProcess([x], [a], [[h]], [f], reshape(σ,1,1), X; time=time, controls=(controls isa Vector ? controls : [controls]))
JumpDiffusionProcess(x::PV, a::Vector{<:APL}, h::Vector{<:APL}, f::APL, σ::APL, X = Fullspace(); time = @polyvar(t)[1], controls=PV{true}[]) =
                     JumpDiffusionProcess([x], a, [[hi] for hi in h], [f], reshape(σ,1,1), X; time=time, controls=(controls isa Vector ? controls : [controls]))

function JumpDiffusionProcess(JP::JumpProcess, DP::DiffusionProcess)
    @assert all(JP.x .== DP.x) "The jump and diffusion process must have the same state"
    return JumpDiffusionProcess(JP.x, JP.a, JP.h, DP.f, DP.σ, intersect(DP.X, JP.X), time = JP.time, controls = JP.controls)
end

mutable struct ControlProcess
    MP::MarkovProcess
    T::Real  # time horizon
    U # set of admissible controls
    Objective # objetive function
    PathChanceConstraints
    TerminalChanceConstraints
    discount_factor
    function ControlProcess(MP::MarkovProcess, T::Real, U, obj, PCs = [], TCs = [], dis_fac = 0)
        return new(MP, T, U, obj, PCs, TCs, dis_fac)
    end
end

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
extended_inf_generator(MP::MarkovProcess, p::Polynomial; scale = 1) = ∂(p,MP.time) + scale*inf_generator(MP, p)
