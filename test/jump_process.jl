using Symbolics, LinearAlgebra, MarkovBounds

function test_JumpProcess(jp, a, h)
    # test propensities
    @test all([isequal(substitute(a[i], jp.poly_vars), jp.a[i]) for i in eachindex(a)])
    # test jumps
    for (k,jump) in enumerate(h)
        @test all([isequal(substitute(jump[i], jp.poly_vars), jp.h[k][i]) for i in eachindex(jump)])
    end
end

@testset "JumpProcess - Tests" begin
    # initialization via Symbolics.jl
    Symbolics.@variables x, y, u, v
    a = [0.2*(x-1.2*v)^2*u; 0.2*y + u]
    h = [[x + 1, y*v], [0.7*x, (x+1.2*y)^2 + u]]
    X = [1.2*(x^2 + 0.2*y)^2 >= 1, (0.2*x - 0.5*y)^2 <= 2.0*x^3]
    jp1 = JumpProcess([x,y], a, h, X, controls = [u,v])
    
    test_JumpProcess(jp1, a, h)
    
    Symbolics.@variables x[1:2], u[1:2]
    a = [0.2*(x[1]-1.2*u[2])^2*u[1]; 0.2*x[2] + u[1]]
    h = [[x[1] + 1, x[2]*u[2]], [0.7*x[1], (x[1]+1.2*x[2])^2 + u[1]]]
    X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
    jp2 = JumpProcess(x, a, h, X, controls = u)

    test_JumpProcess(jp2, a, h)
    
    Symbolics.@variables x[1:2]
    a = [0.2*(x[1]-1.2)^2; 0.2*x[2]]
    h = [[x[1] + 1, x[2]], [0.7*x[1], (x[1]+1.2*x[2])^2]]
    X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
    jp3 = JumpProcess(x, a, h, X)

    test_JumpProcess(jp3, a, h)

    Symbolics.@variables x[1:2], u[1:2]
    x = collect(x)
    u = collect(u)
    a = [0.2*(x[1]-1.2*u[2])^2*u[1]; 0.2*x[2] + u[1]]
    h = [[x[1] + 1, x[2]*u[2]], [0.7*x[1], (x[1]+1.2*x[2])^2 + u[1]]]
    X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
    jp4 = JumpProcess(x, a, h, X, controls = u)

    test_JumpProcess(jp4, a, h)

    # specification via DynamicPolynomials
    @polyvar x[1:2]
    @polyvar u[1:2]
    a = [0.2*(x[1]-1.2*u[2])^2*u[1]; 0.2*x[2] + u[1]]
    h = [[x[1] + 1, x[2]*u[2]], [0.7*x[1], (x[1]+1.2*x[2])^2 + u[1]]]
    X = @set(1.2*(x[1]^2 + 0.2*x[2])^2 >= 1 && (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3)
    jp5 = JumpProcess(x, a, h, X, controls = u)
end