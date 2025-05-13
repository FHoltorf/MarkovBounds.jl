polynomialize_vars(vars::AbstractArray) = Dict(unwrap(v) => _variable(_name(unwrap(v))) for v in vars)
polynomialize_vars(var) = Dict(unwrap(var) => _variable(_name(unwrap(var))))
function polynomialize_vars(x::Arr)  
    x_dp = DynPoly.polyarrayvar(COMMUTATIVE, POLYORDER, x.value.name, x.value.metadata[ArrayShapeCtx]...)
    return Dict(unwrap(x[i]) => x_dp[i] for i in eachindex(x_dp))
end
function polynomialize_vars(vars_args...)
    var_to_poly = Dict()
    for vars in vars_args
        merge!(var_to_poly, polynomialize_vars(vars))
    end
    return var_to_poly
end
#_name(v::Sym) = string(nameof(unwrap(v)))
#_name(v<:BasicSymbolic) = string
#_name(v::Term) = string(v.arguments[1], "[", [string(arg,",") for arg in v.arguments[2:end]]..., "]")
_name(v::Num) = _name(unwrap(v))
function _name(v)
    if istree(v)
        args = arguments(v)
        var_name = string(nameof(args[1]), "[")
        for arg in args[2:end-1]
            var_name *= string(arg, ",")
        end
        var_name *= string(args[end], "]")
        return var_name
    else
        return string(nameof(v))
    end
end

_variable(v::String) = Variable(v, COMMUTATIVE, POLYORDER)

# need to find a better check if expression has been expanded.  
# Here we assume that every expression will be supplied as Num.
# If expression is supplied as NON-EXPANDED SymbolicUtils object, this whole transformation may fail
polynomialize_expr(ex::Num, vars::Dict) = polynomialize_expr(PolyForm(unwrap(expand(ex))), vars)
polynomialize_expr(ex::BasicSymbolic, vars::Dict) = polynomialize_expr(PolyForm(unwrap(expand(ex))), vars)
function polynomialize_expr(pform::PolyForm, vars::Dict)
    w, pvars, psym = pform.p, pform.pvar2sym, pform.sym2term
    Π = [v => (pvars[v] in keys(psym) ? vars[psym[pvars[v]][1]] : vars[pvars[v]]) for v in variables(w)]
    return subs(w, Π...)
end

polynomialize_expr(p::Number, vars::Dict) = polynomial(p)

# polynomialize_expr(ex::Num, vars::Dict) = polynomialize_expr(expand(unwrap(ex)), vars)
# polynomialize_expr(ex::Pow, vars::Dict) = vars[ex.base]^ex.exp
# polynomialize_expr(ex::Mul, vars::Dict) = prod(polynomialize_expr(sub_ex, vars) for sub_ex in arguments(ex))
# polynomialize_expr(ex::Add, vars::Dict) = sum(polynomialize_expr(sub_ex, vars) for sub_ex in arguments(ex))
# polynomialize_expr(ex::Sym, vars::Dict) = vars[ex]
# polynomialize_expr(ex::Term, vars::Dict) = vars[ex]
# polynomialize_expr(ex::Number, ::Dict) = ex

polynomialize_expr(exs::Array, vars::Dict) = POLY[polynomialize_expr(ex, vars) for ex in exs]
polynomialize_expr(exs::SmallVec, vars::Dict) = POLY[polynomialize_expr(ex, vars) for ex in exs]


function polynomialize_set(S::Vector, vars::Dict)
    ineq = POLY[]
    eq = POLY[]
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
    return BasicSemialgebraicSet(algebraic_set(eq), ineq)
end

polynomialize_set(S::Num, vars::Dict) = polynomialize_set([S], vars)

## Process Constructors
# Drift Process
function DriftProcess(x::T, f::Vector{Num}, X = []; iv = [], controls = []) where T <: Union{Arr, Vector{Num}}
    poly_vars = polynomialize_vars(x, controls, iv)
    f = polynomialize_expr(f, poly_vars)
    X_poly = isempty(X) ? FullSpace() : polynomialize_set(X, poly_vars)
    poly_x = [poly_vars[x[i]] for i in eachindex(x)]
    poly_iv = isempty(iv) ? _variable("t") : poly_vars[iv]
    poly_controls = [poly_vars[controls[i]] for i in eachindex(controls)]
    return DriftProcess(poly_x, f, X_poly; poly_vars = poly_vars, iv = poly_iv, controls = poly_controls)
end

# Jump Processes
function JumpProcess(x::T, a::Vector{Num}, h::Vector{Vector{Num}}, X = []; iv = [], controls = []) where T <: Union{Arr, Vector{Num}}
    poly_vars = polynomialize_vars(x, controls, iv)
    a = polynomialize_expr(a, poly_vars)
    h = [polynomialize_expr(hi, poly_vars) for hi in h]
    X_poly = isempty(X) ? FullSpace() : polynomialize_set(X, poly_vars)
    poly_x = [poly_vars[x[i]] for i in eachindex(x)]
    poly_iv = isempty(iv) ? _variable("t") : poly_vars[iv]
    poly_controls = [poly_vars[controls[i]] for i in eachindex(controls)]
    return JumpProcess(poly_x, a, h, X_poly; poly_vars = poly_vars, iv = poly_iv, controls = poly_controls)
end

JumpProcess(x::Num, a::Num, h::Num, X = []; iv = [], controls = []) =
            JumpProcess([x],[a],[[h]],X; iv = iv, controls = controls)
JumpProcess(x::Num, a::Vector{Num}, h::Vector{Num}, X = []; iv = [], controls = []) =
            JumpProcess([x],a,[[hi] for hi in h],X; iv = iv, controls = controls)
JumpProcess(x::Num, a::Vector{Num}, h::Vector{Vector{Num}}, X = []; iv = [], controls = []) =
            JumpProcess([x],a,h,X; iv = iv, controls = controls)

# Diffusion Processes
function DiffusionProcess(x::T, f::Vector{Num}, σ::Matrix{Num}, X = []; iv = [], controls = []) where T <: Union{Vector{Num}, Arr}
    poly_vars = polynomialize_vars(x, controls, iv)
    f = polynomialize_expr(f, poly_vars)
    σ = polynomialize_expr(σ, poly_vars)
    X_poly = isempty(X) ? FullSpace() : polynomialize_set(X, poly_vars)
    poly_x = [poly_vars[x[i]] for i in eachindex(x)]
    poly_iv = isempty(iv) ? _variable("t") : poly_vars[iv]
    poly_controls = [poly_vars[controls[i]] for i in eachindex(controls)]
    return DiffusionProcess(poly_x, f, σ, X_poly, poly_vars = poly_vars, iv = poly_iv, controls = poly_controls)
end

DiffusionProcess(x::Num, f::Num, σ::Num, X = []; iv = [], controls = []) =
                 DiffusionProcess([x], [f], reshape([σ],1,1), X; iv = iv, controls = controls)

# Jump-Diffusion process
function JumpDiffusionProcess(x::T, a::Vector{Num}, h::Vector{Vector{Num}}, f::Vector{Num}, σ::Matrix{Num}, X = []; iv = [], controls = []) where T <: Union{Arr, Vector{Num}}
    poly_vars = polynomialize_vars(x, iv, controls)
    a = polynomialize_expr(a, poly_vars)
    h = [polynomialize_expr(hi, poly_vars) for hi in h]
    f = polynomialize_expr(f, poly_vars)
    σ = polynomialize_expr(σ, poly_vars)
    X_poly = isempty(X) ? FullSpace() : polynomialize_set(X, poly_vars)
    poly_x = [poly_vars[x[i]] for i in eachindex(x)]
    poly_iv = isempty(iv) ? _variable("t") : poly_vars[iv]
    poly_controls = [poly_vars[controls[i]] for i in eachindex(controls)]
    return JumpDiffusionProcess(poly_x, a, h, f, σ, X_poly, poly_vars = poly_vars, iv = poly_iv, controls = poly_controls)
end

JumpDiffusionProcess(x::Num, a::Num, h::Num, f::Num, σ::Num, X = []; iv = [], controls = []) =
                     JumpDiffusionProcess([x], [a], [[h]], [f], reshape([σ],1,1), X; iv = iv, controls = controls)
JumpDiffusionProcess(x::Num, a::Vector{Num}, h::Vector{Num}, f::Num, σ::Num, X = []; iv = [], controls = []) =
                     JumpDiffusionProcess([x], a, [[hi] for hi in h], [f], reshape([σ],1,1), X; iv = iv, controls = controls)
