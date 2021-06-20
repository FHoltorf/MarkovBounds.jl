using Symbolics: Num, expand
using SymbolicUtils: Add, Mul, Pow, Sym, arguments

polynomialize_vars(x, u...) = Dict(v.val => PolyVar{true}(string(v.val.name)) for v in vcat(x,u...))
polynomialize_vars(x) = Dict(v.val => PolyVar{true}(string(v.val.name)) for v in x)

# need to find a better check if expression has been expanded. Here I assume that every expression will be supplied as Num.
# If expression is supplied as NON-EXPANDED SymbolicUtils object, this whole transformation may fail
# For now good enough.
polynomialize_expr(ex::Num, vars::Dict) = polynomialize_expr(expand(ex).val, vars)
polynomialize_expr(ex::Pow, vars::Dict) = vars[ex.base]^ex.exp
polynomialize_expr(ex::Mul, vars::Dict) = prod(polynomialize_expr(sub_ex, vars) for sub_ex in arguments(ex))
polynomialize_expr(ex::Add, vars::Dict) = sum(polynomialize_expr(sub_ex, vars) for sub_ex in arguments(ex))
polynomialize_expr(ex::Sym, vars::Dict) = vars[ex]
polynomialize_expr(ex::Number, vars::Dict) = ex

function polynomialize_expr(exs::Array, vars::Dict)
    poly_exs = Array{Polynomial{true, Float64}}(undef, size(exs))
    for (i, ex) in enumerate(exs)
        poly_exs[i] = polynomialize_expr(ex, vars)
    end
    return poly_exs
end

function polynomialize_set(S::Vector, vars::Dict)
    ineq = Polynomial{true,Float64}[]
    eq = Polynomial{true,Float64}[]
    for entry in S
        polys = polynomialize_expr(entry.val.arguments, vars)
        if entry.val.f == >=
            push!(ineq, polys[1] - polys[2])
        elseif entry.val.f == <=
            push!(ineq, polys[2] - polys[1])
        else
            push!(eq, polys[1] - polys[2])
        end
    end
    return BasicSemialgebraicSet(algebraicset(eq), ineq)
end

polynomialize_set(S::Num, vars::Dict) = polynomialize_set([S], vars)


## Process Constructors
# Jump Processes
function JumpProcess(x::Vector{Num}, a::Vector{Num}, h::Vector{Vector{Num}}, X = []; time = [], controls = [])
    poly_vars = polynomialize_vars(x, controls, time)
    a = polynomialize_expr(a, poly_vars)
    h = [polynomialize_expr(hi, poly_vars) for hi in h]
    X_poly = isempty(X) ? FullSpace() : polynomialize_set(X, poly_vars)
    poly_x = [poly_vars[state] for state in x]
    poly_time = isempty(time) ? @polyvar(t)[1] : poly_vars[time]
    poly_controls = isempty(controls) ? PV{true}[] : [poly_vars[u] for u in controls]
    return JumpProcess(poly_x, a, h, X_poly; poly_vars = poly_vars, time = poly_time, controls = poly_controls)
end

JumpProcess(x::Num, a::Num, h::Num, X = []; time = [], controls = []) =
            JumpProcess([x],[a],[[h]],X; time = time, controls = controls)
JumpProcess(x::Num, a::Vector{Num}, h::Vector{Num}, X = []; time = [], controls = []) =
            JumpProcess([x],a,[[hi] for hi in h],X; time = time, controls = controls)
JumpProcess(x::Num, a::Vector{Num}, h::Vector{Vector{Num}}, X = []; time = [], controls = []) =
            JumpProcess([x],a,h,X; time = time, controls = controls)

# Diffusion Processes
function DiffusionProcess(x::Vector{Num}, f::Vector{Num}, σ::Matrix{Num}, X = []; time = [], controls = [])
    poly_vars = polynomialize_vars(x, controls, time)
    f = polynomialize_expr(f, poly_vars)
    σ = polynomialize_expr(σ, poly_vars)
    X_poly = isempty(X) ? FullSpace() : polynomialize_set(X, poly_vars)
    poly_x = [poly_vars[state] for state in x]
    poly_time = isempty(time) ? @polyvar(t)[1] : poly_vars[time]
    poly_controls = isempty(controls) ? PV{true}[] : [poly_vars[u] for u in controls]
    return DiffusionProcess(poly_x, f, σ, X_poly, poly_vars = poly_vars, time = poly_time, controls = poly_controls)
end

DiffusionProcess(x::Num, f::Num, σ::Num, X = []; time = [], controls = []) =
                 DiffusionProcess([x], [f], reshape([σ],1,1), X; time = time, controls = controls)

# Jump-Diffusion process
function JumpDiffusionProcess(x::Vector{Num}, a::Vector{Num}, h::Vector{Vector{Num}}, f::Vector{Num}, σ::Matrix{Num}, X = []; time = [], controls = [])
    poly_vars = polynomialize_vars(x, time, controls)
    a = polynomialize_expr(a, poly_vars)
    h = [polynomialize_expr(hi, poly_vars) for hi in h]
    f = polynomialize_expr(f, poly_vars)
    σ = polynomialize_expr(σ, poly_vars)
    X_poly = isempty(X) ? FullSpace() : polynomialize_set(X, poly_vars)
    poly_x = [poly_vars[state] for state in x]
    poly_time = isempty(time) ? @polyvar(t)[1] : poly_vars[time]
    poly_controls = isempty(controls) ? PV{true}[] : [poly_vars[u] for u in controls]
    return JumpDiffusionProcess(poly_x, a, h, f, σ, X_poly, poly_vars = poly_vars, time = poly_time, controls = poly_controls)
end

JumpDiffusionProcess(x::Num, a::Num, h::Num, f::Num, σ::Num, X = []; time = [], controls = []) =
                     JumpDiffusionProcess([x], [a], [[h]], [f], reshape([σ],1,1), X; time = time, controls = controls)
JumpDiffusionProcess(x::Num, a::Vector{Num}, h::Vector{Num}, f::Num, σ::Num, X = []; time = [], controls = []) =
                     JumpDiffusionProcess([x], a, [[hi] for hi in h], [f], reshape([σ],1,1), X; time = time, controls = controls)
