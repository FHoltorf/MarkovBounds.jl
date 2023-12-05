function test_DriftProcess(dp, f)
    # test drift
    @test all([isequal(Symbolics.substitute(f[i], dp.poly_vars), dp.f[i]) for i in eachindex(f)])
end

# initialization via Symbolics.jl
@testset "DriftProcess - Constructor Tests" begin
    Symbolics.@variables x, y, u, v
    f = [1.5*x*(x-1)^2*u-v^2; 1.2*y+(u+0.2)^2]
    X = [1.2*(x^2 + 0.2*y)^2 >= 1, (0.2*x - 0.5*y)^2 <= 2.0*x^3]
    dp1 = DriftProcess([x,y], f, X, controls = [u,v])

    test_DriftProcess(dp1, f)

    Symbolics.@variables x[1:2], u[1:2]
    f = [1.5*x[1]*(x[1]-1)^2*u[1]-u[2]^2; 1.2*x[2]+(u[1]+0.2)^2]
    X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
    dp2 = DriftProcess(x, f, X, controls = u)

    test_DriftProcess(dp2, f)
    
    Symbolics.@variables x[1:2]
    f = [1.5*x[1]*(x[1]-1)^2; 1.2*x[2]]
    X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
    dp3 = DriftProcess(x, f, X)

    test_DriftProcess(dp3, f)
    
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

@testset "DriftProcess - Optimization test" begin
    @polyvar(x[1:2]) # state variables
    @polyvar(u) # control variables
    @polyvar(t) # time variable
    X = @set(x[1] >= 0 && x[2] >= 0) # state space

    γ = [1, 2, 1, 2, 0.25*0.1] # model parameters

    f = [γ[1] * x[1] - γ[2] * x[1] * x[2] ;
        γ[4] * x[1] * x[2] - γ[3] * x[2] - x[2]*u] # drift coefficient

    lv = DriftProcess(x, f, X, iv = t, controls = [u]) 

    U = @set(u >= 0 && u <= 1) # set of admissible controls|
    stagecost = (x[1]-0.75)^2 + (x[2] - 0.5)^2/10 + (u - 0.5)^2/10
    obj = Lagrange(stagecost) # Lagrange type objective
    T = 10.0 # control horizon
    lv_control = ControlProcess(lv, T, U, obj)

    order = 2
    x0 = [1.0, 0.25]
    μ0 = Dict(x[1]^i*x[2]^j => x0[1]^i*x0[2]^j for i in 0:order+1, j in 0:order+1) # moments of initial distribution
    trange = range(0, T, length = 11) # discretization of time horizon
    b = optimal_control(lv_control, μ0, order, trange, solver)
    @test b.value >= 0.0
end