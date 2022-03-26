using MarkovBounds, COSMO, Catalyst
using Test

opt = COSMO
@testset "stationary bounds" begin
    # Birth Death Process
    BD = @reaction_network begin
        1.0, ∅ --> A
        0.01, 2A --> A
    end
    A = species(BD)[1]
    A0 = Dict(A => 2.0)
    A_scale = Dict(A => 10.0)
    lb, ub = stationary_mean(BD, A0, A, 4, opt.Optimizer, A_scale)
    @test lb.value <= ub.value

    ub = stationary_variance(BD, A0, A, 2, opt.Optimizer, A_scale)
    @test ub.value >= 0

    # Michaelis-Menten System
    MM = @reaction_network begin
        10.0, S + E --> SE
        0.2, SE --> S + E
        4.0, SE --> P + E
        1.3, P --> S
    end
    S, E, SE, P = species(MM)
    x0 = Dict(S => 10, E => 5, SE => 0, P => 0)
    x_scale = Dict(S => 10, E => 10, SE => 10, P => 10)
    lb = stationary_covariance_ellipsoid(MM, x0, [S, P], 2, opt.Optimizer, x_scale)
    @test lb.value >= 0
end

@testset "transient bounds" begin
    BD = @reaction_network begin
        1.0, ∅ --> A
        0.01, 2A --> A
    end
    A = species(BD)[1]
    A0 = Dict(A => 2.0)
    A_scale = Dict(A => 10.0)
    trange = 0:0.1:1.0
    lb, ub = transient_mean(BD, A0, A, 2, trange, opt.Optimizer, A_scale)
    @test lb.value <= ub.value

    ub = transient_variance(BD, A0, A, 2, trange, opt.Optimizer, A_scale)
    @test ub.value >= 0

    # Michaelis-Menten System
    MM = @reaction_network begin
        10.0, S + E --> SE
        0.2, SE --> S + E
        4.0, SE --> P + E
        1.3, P --> S
    end
    S, E, SE, P = species(MM)
    x0 = Dict(S => 10, E => 5, SE => 0, P => 0)
    x_scale = Dict(S => 10, E => 10, SE => 10, P => 10)
    ub = transient_covariance_ellipsoid(MM, x0, [S, P], 2, trange, opt.Optimizer, x_scale)
    @test ub.value >= 0
end

@testset "Optimal Control" begin
    @polyvar(x[1:2])
    @polyvar(u[1:1])
    @polyvar(t)
    γ = [1, 2, 1, 2, 0.25*0.1]

    Tf = 10.0
    x_target = [0.75, 0.5]
    u_target = 0.5

    f = [γ[1] * x[1] - γ[2] * x[1] * x[2] ;
         γ[4] * x[1] * x[2] - γ[3] * x[2] - x[2]*u[1]]

    g = [γ[5]*x[1]; 0]
    σ = polynomial.(g*g')

    lagrange = (x[1] - x_target[1])^2 + (x[2] - x_target[2])^2/10 + (u[1] - u_target)^2/10

    trange = 1.0:1.0:Tf
    MP = DiffusionProcess(x, f, σ, @set(x[1] >= 0 && x[2] >= 0), iv = t, controls = u)
    CP = ControlProcess(MP, Tf, u, t, @set(u[1] >= 0 && u[1] <= 1), LagrangeMayer(lagrange, 0*x[1]))
    @test CP.MP.controls == u
end
