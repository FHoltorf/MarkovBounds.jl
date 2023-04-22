using Symbolics, LinearAlgebra, MarkovBounds

function test_DriftProcess(dp, f)
    # test drift
    @test all([isequal(substitute(f[i], dp.poly_vars), dp.f[i]) for i in eachindex(f)])
end

# initialization via Symbolics.jl
@testset "DriftProcess - tests" begin
    Symbolics.@variables x, y, u, v
    f = [1.5*x*(x-1)^2*u-v^2; 1.2*y+(u+0.2)^2]
    X = [1.2*(x^2 + 0.2*y)^2 >= 1, (0.2*x - 0.5*y)^2 <= 2.0*x^3]
    dp1 = DriftProcess([x,y], f, X, controls = [u,v])

    test_DriftProcess(dp1, f)

    Symbolics.@variables x[1:2], u[1:2]
    f = [1.5*x[1]*(x[1]-1)^2*u[1]-u[2]^2; 1.2*x[2]+(u[1]+0.2)^2]
    σ = Matrix(Diagonal([0.2*x[1]^2; 0.2*x[2]]))
    X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
    dp2 = DriftProcess(x, f, X, controls = u)

    test_DriftProcess(dp2, f)
    
    Symbolics.@variables x[1:2]
    f = [1.5*x[1]*(x[1]-1)^2; 1.2*x[2]]
    σ = Matrix(Diagonal([0.2*x[1]^2; 0.2*x[2]]))
    X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
    dp3 = DriftProcess(x, f, X)

    test_DriftProcess(dp3, f, σ)
    
    Symbolics.@variables x[1:2], u[1:2]
    x = collect(x)
    u = collect(u)
    f = [1.5*x[1]*(x[1]-1)^2*u[1]-u[2]^2; 1.2*x[2]+(u[1]+0.2)^2]
    X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
    dp4 = DriftProcess(x, f, X, controls = u)

    test_DriftProcess(dp4, f)
    
    # specification via DynamicPolynomials
    @polyvar x[1:2]
    @polyvar u[1:2]
    f = [1.5*x[1]*(x[1]-1)^2*u[1]-u[2]^2; 1.2*x[2]+(u[1]+0.2)^2]
    X = @set(0.2*(x[1]^2 + 0.2*x[2])^2 - 1 >= 0 && (0.2*x[1] - 0.5*x[2])^2 - 2.0*x[1]^3 >= 0)
    dp4 = DriftProcess(x, f, X, controls = u)

    test_DriftProcess(dp4, f)
end