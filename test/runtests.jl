using MarkovBoundsSOS, MosekTools
using Test

@testset "stationary bounds" begin
    # Birth Death Process
    BD = @reaction_network begin
        1.0, ∅ --> A
        0.01, 2A --> A
    end
    A = species(BD)[1]
    A0 = Dict(A => 2.0)
    A_scale = Dict(A => 10.0)
    val1, stat, time = stationary_mean(BD, A0, A, 2, Mosek.Optimizer, A_scale)
    val2, stat, time = stationary_mean(BD, A0, A, 4, Mosek.Optimizer, A_scale)

    @test val1[1] <= val1[2] #mean bounds - validity
    @test val2[1] >= val1[1] && val2[2] <= val1[2] #mean bounds - tightening

    val1, stat, time = stationary_variance(BD, A0, A, 2, Mosek.Optimizer, A_scale)
    val2, stat, time = stationary_variance(BD, A0, A, 4, Mosek.Optimizer, A_scale)

    @test val1 >= 0 #variance bound - validity
    @test val2 <= val1 #variance bound - tightening

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
    val1, stat, time = stationary_covariance_ellipsoid(MM, x0, [S, P], 2, Mosek.Optimizer, x_scale)
    val2, stat, time = stationary_covariance_ellipsoid(MM, x0, [S, P], 4, Mosek.Optimizer, x_scale)
    @test val1 >= 0 # covariance ellipsoid - validity
    @test val2 <= val1 #covariance ellipsoid - tightening
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
    val1, stat, time = transient_mean(BD, A0, A, 2, trange, Mosek.Optimizer, A_scale)
    val2, stat, time = transient_mean(BD, A0, A, 4, trange, Mosek.Optimizer, A_scale)
    @test val1[1] <= val1[2] #mean bounds - validity
    @test val2[1] >= val1[1] && val2[2] <= val1[2] #mean bounds - tightening

    val1, stat, time = transient_variance(BD, A0, A, 2, trange, Mosek.Optimizer, A_scale)
    val2, stat, time = transient_variance(BD, A0, A, 4, trange, Mosek.Optimizer, A_scale)

    @test val1 >= 0 #variance bound - validity
    @test val2 <= val1 #variance bound - tightening

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
    val1, stat, time = transient_covariance_ellipsoid(MM, x0, [S, P], 2, trange, Mosek.Optimizer, x_scale)
    val2, stat, time = transient_covariance_ellipsoid(MM, x0, [S, P], 4, trange, Mosek.Optimizer, x_scale)
    @test val1 >= 0 # covariance ellipsoid - validity
    @test val2 <= val1 #covariance ellipsoid - tightening
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

    trange = 0:1.0:Tf
    MP = DiffusionProcess(x, f, σ, @set(x[1] >= 0 && x[2] >= 0))
    CP = ControlProcess(MP, Tf, u, t, @set(u[1] >= 0 && u[1] <= 1), LagrangeMayer(lagrange, 0*x[1]))

    x0 = [1.0; 0.25]
    μ0 = Dict(x[1]^i*x[2]^j => x0[1]^i*x0[2]^j for i in 0:10, j in 0:10)
    res1 = optimal_control(CP, μ0, 2, trange, Mosek.Optimizer)
    res2 = optimal_control(CP, μ0, 4, trange, Mosek.Optimizer)
    @test res1 <= res2
end
