polynomialize_vars(vars::AbstractArray) = Dict(v.val => PolyVar{true}(name(v.val)) for v in vars)
polynomialize_vars(var) = Dict(var.val => PolyVar{true}(name(v.val)))
function polynomialize_vars(x::Arr)  
    x_dp = polyarrayvar(PolyVar{true}, x.value.name, x.value.metadata[ArrayShapeCtx]...)
    return Dict(x[i].val => x_dp[i] for i in 1:length(x_dp))
end
function polynomialize_vars(vars_args...)
    var_to_poly = Dict()
    for vars in vars_args
        merge!(var_to_poly, polynomialize_vars(vars))
    end
    return var_to_poly
end
name(v::Sym) = string(v.name)
name(v::Term) = string(v.arguments[1], "[", [string(arg,",") for arg in v.arguments[2:end]]..., "]")

# need to find a better check if expression has been expanded. 
# Here we assume that every expression will be supplied as Num.
# If expression is supplied as NON-EXPANDED SymbolicUtils object, this whole transformation may fail
polynomialize_expr(ex::Num, vars::Dict) = polynomialize_expr(expand(ex).val, vars)
polynomialize_expr(ex::Pow, vars::Dict) = vars[ex.base]^ex.exp
polynomialize_expr(ex::Mul, vars::Dict) = prod(polynomialize_expr(sub_ex, vars) for sub_ex in arguments(ex))
polynomialize_expr(ex::Add, vars::Dict) = sum(polynomialize_expr(sub_ex, vars) for sub_ex in arguments(ex))
polynomialize_expr(ex::Sym, vars::Dict) = vars[ex]
polynomialize_expr(ex::Term, vars::Dict) = vars[ex]
polynomialize_expr(ex::Number, ::Dict) = ex

function polynomialize_expr(exs::Array, vars::Dict)
    poly_exs = Array{Polynomial{true, Float64}}(undef, size(exs))
    for (i, ex) in enumerate(exs)
        poly_exs[i] = polynomialize_expr(expand(ex), vars)
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
function JumpProcess(x::T, a::Vector{Num}, h::Vector{Vector{Num}}, X = []; iv = [], controls = []) where T <: Union{Arr, Vector{Num}}
    poly_vars = polynomialize_vars(x, controls, iv)
    a = polynomialize_expr(a, poly_vars)
    h = [polynomialize_expr(hi, poly_vars) for hi in h]
    X_poly = isempty(X) ? FullSpace() : polynomialize_set(X, poly_vars)
    poly_x = PV{true}[poly_vars[x[i]] for i in eachindex(x)]
    poly_iv = isempty(iv) ? PV{true}("t") : poly_vars[iv]
    poly_controls = PV{true}[poly_vars[controls[i]] for i in eachindex(controls)]
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
    poly_iv = isempty(iv) ? PV{true}("t") : poly_vars[iv]
    poly_controls = PV{true}[poly_vars[controls[i]] for i in eachindex(controls)]
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
    poly_iv = isempty(iv) ? PV{true}("t") : poly_vars[iv]
    poly_controls = [poly_vars[u[i]] for i in eachindex(controls)]
    return JumpDiffusionProcess(poly_x, a, h, f, σ, X_poly, poly_vars = poly_vars, iv = poly_iv, controls = poly_controls)
end

JumpDiffusionProcess(x::Num, a::Num, h::Num, f::Num, σ::Num, X = []; iv = [], controls = []) =
                     JumpDiffusionProcess([x], [a], [[h]], [f], reshape([σ],1,1), X; iv = iv, controls = controls)
JumpDiffusionProcess(x::Num, a::Vector{Num}, h::Vector{Num}, f::Num, σ::Num, X = []; iv = [], controls = []) =
                     JumpDiffusionProcess([x], a, [[hi] for hi in h], [f], reshape([σ],1,1), X; iv = iv, controls = controls)
