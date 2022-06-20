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