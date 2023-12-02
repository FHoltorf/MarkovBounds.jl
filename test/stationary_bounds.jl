using DynamicPolynomials, Catalyst, COSMO
# specify solver
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
lotka_volterra = DiffusionProcess(x, f, σ, X)

# define simple birth-death process
birth_death = @reaction_network begin
                1.0, ∅ --> A
                0.01, 2A --> A
               end
A = species(birth_death)[1]
A0 = Dict(A => 2.0)
birth_death_jump = ReactionProcess(birth_death, A0)
birth_death_langevin = LangevinProcess(birth_death, A0)

@testset "Stationary bounds for reaction processes" begin
    # stationary_polynomial problem 
    val = 15.1421
    bound = stationary_polynomial(birth_death_jump, Num(A+1), 2, solver)
    @test abs(bound.value - val) < 1e-3
    bound = stationary_polynomial(birth_death_langevin, Num(A+1), 2, solver)
    @test abs(bound.value - val) < 1e-3
    
    # stationary_mean problem
    val -= 1
    lb, ub = stationary_mean(birth_death_jump, A, 2, solver)
    @test abs(lb.value - val) < 1e-3
    lb, ub = stationary_mean(birth_death_langevin, A, 2, solver)
    @test abs(lb.value - val) < 1e-3
    lb, ub = stationary_mean(birth_death, A0, A, 2, solver)
    @test abs(lb.value - val) < 1e-3
    lb, ub = stationary_mean(birth_death, A, 2, solver)
    @test abs(lb.value - val) < 1e-3
    
    # stationary_variance problem 
    val = 14.14227
    ub = stationary_variance(birth_death_jump, A, 2, solver)
    @test abs(ub.value - val) < 1e-3
    ub = stationary_variance(birth_death_langevin, A, 2, solver)
    @test abs(ub.value - val) < 1e-3
    ub =stationary_variance(birth_death, A, 2, solver)    
    @test abs(ub.value - val) < 1e-3
end
    
@testset "stationary_polynomial tests" begin
    lb = stationary_polynomial(lotka_volterra, x[1] + x[2], 2, solver)
    @test (lb.value - 0.0) <= 1e-3
end

@testset "stationary_covariance_ellipsoid tests" begin
    lb = stationary_covariance_ellipsoid(lotka_volterra, [x[1], x[2]], 2, solver)
    @test (lb.value - 0.0) <= 1e-3    
end

@testset "stationary_variance tests" begin
    ub = stationary_variance(lotka_volterra, x[1]+x[2], 2, solver)
    @test ub.value >= 0
end

@testset "stationary_probability_mass tests" begin
    S = @set(x[1] <= 0.01 && x[2] <= 0.01)
    lb, ub = MarkovBounds.stationary_probability_mass(lotka_volterra, S, 2, solver)
    println(lb.value)
    println(ub.value)
    @test lb.value >= -1e-6
    @test ub.value <= 1 + 1e6
end

#TBD
#=
@testset "approximate_stationary_measure" begin

end

@testset "max_entropy_measure" begin

end
=#
