using Symbolics: Num, expand
using SymbolicUtils: Mul, Pow, Sym, arguments

polynomialize_vars(x, u) = Dict(v.val => PolyVar{true}(string(v.val.name)) for v in vcat(x,u))
polynomialize_vars(x) = Dict(v.val => PolyVar{true}(string(v.val.name)) for v in x)

# one can probably do this better via multiple dispatch but for now it does what I want
function polynomialize_expr(ex, vars)
    ex = SymbolicUtils.expand(ex)
    if ex isa SymbolicUtils.Mul
        poly = parse_prod(SymbolicUtils.arguments(ex), vars)
    elseif ex isa SymbolicUtils.Add
        poly = parse_sum(SymbolicUtils.arguments(ex), vars)
    else
        @error "Expression $ex cannot be transformed into polynomial type. Potentially bug"
    end
    return poly
end

function parse_sum(args,  vars)
    poly = 0
    for tm in args
        if tm isa Number
            poly += tm
        elseif tm isa SymbolicUtils.Mul
            poly += tm.coeff*prod(vars[var]^tm.dict[var] for var in keys(tm.dict))
        elseif tm isa SymbolicUtils.Pow
            @assert tm.exp isa Int "only polynomials are supported"
            poly += vars[tm.base]^tm.exp
        elseif tm isa SymbolicUtils.Sym
            poly += vars[tm]
        else
            @error "unsupported component in function"
        end
    end
    return poly
end

function parse_prod(args, vars)
    poly = 1
    for tm in args
        if tm isa Number
            poly *= tm
        elseif tm isa Mul
            poly *= tm.coeff*prod(vars[var]^tm.dict[var] for var in keys(tm.dict))
        elseif tm isa Pow
            @assert tm.exp isa Int "only polynomials are supported"
            poly *= vars[tm.base]^tm.exp
        elseif tm isa Sym
            poly *= vars[tm]
        else
            @error "unsupported component in function"
        end
    end
    return poly
end

polynomialize_expr(ex::Sym, vars::Dict) = vars[ex]
polynomialize_expr(ex::Num, vars::Dict) = polynomialize_expr(ex.val, vars)
polynomialize_expr(ex::Number, ::Dict) = ex

function polynomialize_expr(exs::T, vars::Dict) where T <: Union{Vector, Matrix}
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
function JumpProcess(x::Vector{Num}, a::Vector{Num}, h::Vector{Vector{Num}}, X = []; u = [])
    poly_vars = polynomialize_vars(x, u)
    a = polynomialize_expr(a, poly_vars)
    h = [polynomialize_expr(hi, poly_vars) for hi in h]
    X_poly = isempty(X) ? FullSpace() : polynomialize_set(X, poly_vars)
    states = [poly_vars[state] for state in x]
    return JumpProcess(states, a, h, X_poly, poly_vars = poly_vars)
end

JumpProcess(x::Num, a::Num, h::Num, X = []; u = []) = JumpProcess([x],[a],[[h]],X; u = u)
JumpProcess(x::Num, a::Vector{Num}, h::Vector{Num}, X = []; u = []) = JumpProcess([x],a,[[hi] for hi in h],X; u = u)
JumpProcess(x::Num, a::Vector{Num}, h::Vector{Vector{Num}}, X = []; u = []) = JumpProcess([x],a,h,X; u = u)

# Diffusion Processes
function DiffusionProcess(x::Vector{Num}, f::Vector{Num}, σ::Matrix{Num}, X = []; u = [])
    poly_vars = polynomialize_vars(x, u)
    f = polynomialize_expr(f, poly_vars)
    σ = polynomialize_expr(σ, poly_vars)
    X_poly = isempty(X) ? FullSpace() : polynomialize_set(X, poly_vars)
    states = [poly_vars[state] for state in x]
    return DiffusionProcess(states, f, σ, X_poly, poly_vars = poly_vars)
end

DiffusionProcess(x::Num, f::Num, σ::Num, X = []; u = []) = DiffusionProcess([x], [f], reshape([σ],1,1), X; u = u)

# Jump-Diffusion process
function JumpDiffusionProcess(x::Vector{Num}, a::Vector{Num}, h::Vector{Vector{Num}}, f::Vector{Num}, σ::Matrix{Num}, X = []; u = [])
    poly_vars = polynomialize_vars(x, u)
    a = polynomialize_expr(a, poly_vars)
    h = [polynomialize_expr(hi, poly_vars) for hi in h]
    f = polynomialize_expr(f, poly_vars)
    σ = polynomialize_expr(σ, poly_vars)
    X_poly = isempty(X) ? FullSpace() : polynomialize_set(X, poly_vars)
    states = [poly_vars[state] for state in x]
    return JumpDiffusionProcess(states, a, h, f, σ, X_poly, poly_vars = poly_vars)
end

JumpDiffusionProcess(x::Num, a::Num, h::Num, f::APL, σ::APL, X = []; u = []) = JumpDiffusionProcess([x], [a], [[h]], [f], reshape(σ,1,1), X; u = u)
JumpDiffusionProcess(x::Num, a::Vector{Num}, h::Vector{Num}, f::Num, σ::Num, X = []; u = []) = JumpDiffusionProcess([x], a, [[hi] for hi in h], [f], reshape(σ,1,1), X; u = u)
