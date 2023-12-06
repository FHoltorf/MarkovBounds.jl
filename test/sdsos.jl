#=
# test not interesting
birth_death = @reaction_network begin
    1.0, ∅ --> A
    0.01, 2A --> A
end
A = species(birth_death)[1]
A0 = Dict(A => 2.0)
A_scale = Dict(A => 10.0)

rp, A0_scaled = reaction_process_setup(birth_death, A0; scales = A_scale)
for d in [2,4,6,8]
    model = SOSModel(solver)
    PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly,SOSCone)
    monos = monomials(rp.JumpProcess.x, 0:d)
    @variable(model, w, Poly(monos))
    @variable(model, γ)
    @constraint(model, MarkovBounds.inf_generator(rp.JumpProcess, w) + rp.JumpProcess.x[1] >= γ, 
                       domain = rp.JumpProcess.X, maxdegree=d+10)
    @objective(model, Max, γ)
    optimize!(model)
    println("sos: $(objective_value(model))")
end
=#


@testset "(S)DSOS hierarchies" begin
    # bounding polynomials
    @polyvar(x[1:2]) # state variables
    u = 0.5 # constant hunting 
    X = @set(x[1] >= 0 && x[2] >= 0) # state space
    γ = [1, 2, 1, 2, 0.25*0.1] # model parameters
    f = [γ[1] * x[1] - γ[2] * x[1] * x[2] ;
        γ[4] * x[1] * x[2] - γ[3] * x[2] - x[2]*u] 
    g = [γ[5]*x[1]; 0] # diffusion coefficient
    σ = polynomial.(g*g') 
    lotka_volterra = DiffusionProcess(x, f, σ, X)
    l = - (x[1]^2 + x[2]^2)

    for cone in [SDSOSCone, SOSCone]
        for d in [4, 6]
            model = SOSModel(solver)
            PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, cone)
            monos = monomials(x, 0:d)
            @variable(model, w, Poly(monos))
            @variable(model, γ)
            @constraint(model, MarkovBounds.inf_generator(lotka_volterra, w) + l >= γ, domain = lotka_volterra.X)
            @objective(model, Max, γ)
            optimize!(model)
            b = stationary_polynomial(lotka_volterra, l, d, solver; inner_approx=cone)
            @test abs(b.value - objective_value(model))/b.value <= 0.01
        end
    end

    @polyvar u 
    X = @set(x[1] >= 0 && x[2] >= 0) # state space
    U = @set(u >= 0 && u <= 1)
    γ = [1, 2, 1, 2, 0.25*0.1] # model parameters
    f = [γ[1] * x[1] - γ[2] * x[1] * x[2] ;
        γ[4] * x[1] * x[2] - γ[3] * x[2] - x[2]*u] 
    g = [γ[5]*x[1]; 0] # diffusion coefficient
    σ = polynomial.(g*g') 
    lotka_volterra = DiffusionProcess(x, f, σ, X, controls = [u])

    U = @set(u >= 0 && u <= 1) # set of admissible controls|
    stagecost = (x[1]-0.75)^2 + (x[2] - 0.5)^2/10 + (u - 0.5)^2/10
    obj = Lagrange(stagecost) # Lagrange type objective
    T = 10 # control horizon

    inf_horizon_control = ControlProcess(lotka_volterra, Inf, U, obj, [], [], 0.1)
    x0 = [1.0, 0.25]
    for d in [2,4,6,8]
        μ0 = init_moments(x, x0, d)
        sos = optimal_control(inf_horizon_control, μ0, d, 0:5.0:T, solver)
        sdsos = optimal_control(inf_horizon_control, μ0, d, 0:5.0:T, solver, inner_approx=SDSOSCone)
        dsos = optimal_control(inf_horizon_control, μ0, d, 0:5.0:T, solver, inner_approx=DSOSCone) 
        @test sos.value >= sdsos.value >= dsos.value
        println("sos: $(sos.value), sdsos: $(sdsos.value), dsos: $(dsos.value)")
    end
end