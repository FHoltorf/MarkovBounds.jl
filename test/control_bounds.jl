using DynamicPolynomials, Catalyst, MarkovBounds, COSMO
# define solver
solver = COSMO.Optimizer #Hypatia.Optimizer

# define diffusion process -> Lotka-Volterra Predator-Prey model
@polyvar(x[1:2]) # state variables
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

@testset "control problem tests" begin
    fin_horizon_control = ControlProcess(lotka_volterra, T, U, obj)
    inf_horizon_control = ControlProcess(lotka_volterra, Inf, U, obj, [], [], 0.1)

    x0 = [1.0, 0.25]
    μ0 = init_moments(x, x0, 4)
    inf_horizon = optimal_control(inf_horizon_control, μ0, 2, 0:5.0:T, solver)
    @test inf_horizon.value > 0

    fin_horizon = optimal_control(fin_horizon_control, μ0, 2, 0:5.0:T, solver)
    @test fin_horizon.value > 0
end