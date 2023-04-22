using DynamicPolynomials, Catalyst, COSMO
# define solver
solver = COSMO.Optimizer

# define diffusion process -> Lotka-Volterra Predator-Prey model
@polyvar(x[1:2]) # state variables
u = 0.5 # constant hunting 
X = @set(x[1] >= 0 && x[2] >= 0) # state space
γ = [1, 2, 1, 2, 0.25*0.1] # model parameters
f = [γ[1] * x[1] - γ[2] * x[1] * x[2] ;
     γ[4] * x[1] * x[2] - γ[3] * x[2] - x[2]*u] 
g = [γ[5]*x[1]; 0] # diffusion coefficient
σ = polynomial.(g*g') 
x0 = [0.7, 0.5]
lotka_volterra = DiffusionProcess(x, f, σ, X)