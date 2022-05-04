using MarkovBounds, Catalyst
using Test

@testset "Reaction Process - Tests" begin
    # Birth Death Process
    BD = @reaction_network begin
        1.0, ∅ --> A
        0.01, 2A --> A
    end
    A = species(BD)[1]
    A0 = Dict(A => 2.0)
    A_scale = Dict(A => 10.0)
    
    rp = ReactionProcess(BD)
    x = rp.JumpProcess.x
    @test all(rp.JumpProcess.a .== [1.0, 0.01*x[1]*(x[1]-1)/2])
    @test all(rp.JumpProcess.h .== [[x[1]+1], [x[1]-1]])
    @test all(rp.JumpProcess.X.p == [x[1]])
    
    rp, A0_scaled = reaction_process_setup(BD, A0; scales = A_scale)
    x = rp.JumpProcess.x
    @test A0_scaled[1] == A0[A]/A_scale[A]
    @test all(rp.JumpProcess.a .== [1.0, 0.01*x[1]*A_scale[A]*(A_scale[A]*x[1]-1)/2])
    @test all(rp.JumpProcess.h .== [[x[1]+1/A_scale[A]], [x[1]-1/A_scale[A]]])
    @test all(rp.JumpProcess.X.p .== [A_scale[A]*x[1]])
    
    # Michaelis-Menten System
    MM = @reaction_network begin
        10.0, S + E --> SE
        0.2, SE --> S + E
        4.0, SE --> P + E
        1.3, P --> S
    end
    S, E, SE, P = species(MM)
    x0 = Dict(S => 10, E => 5, SE => 0, P => 0)
    x_scale = Dict(S => 5, E => 6, SE => 7, P => 8)

    rp = ReactionProcess(MM)
    x = rp.JumpProcess.x
    @test all(rp.JumpProcess.a .== [10.0*rp.species_to_state[S]*rp.species_to_state[E],
                                    0.2*rp.species_to_state[SE],
                                    4.0*rp.species_to_state[SE],
                                    1.3*rp.species_to_state[P]])
    @test all(rp.JumpProcess.h .== [x .+ [-1, -1, 1, 0], 
                                    x .+ [1, 1, -1, 0],
                                    x .+ [0, 1, -1, 1],
                                    x .+ [1, 0, 0, -1]])
    @test all(rp.JumpProcess.X.p .== [x[1], x[2], x[3], x[4]]) 

    rp, x0_scaled = reaction_process_setup(MM, x0; scales = x_scale)
    x = rp.JumpProcess.x
    @test all(x0_scaled .== [0.0,0.0])
    @test all(rp.JumpProcess.a .== [10.0*rp.species_to_state[S]*rp.species_to_state[E],
                                    0.2*rp.species_to_state[SE],
                                    4.0*rp.species_to_state[SE],
                                    1.3*rp.species_to_state[P]])
    @test all(rp.JumpProcess.h .== [x .+ [1/x_scale[SE], 0],
                                    x .+ [-1/x_scale[SE], 0],
                                    x .+ [-1/x_scale[SE], 1/x_scale[P]],
                                    x .+ [0, -1/x_scale[P]]])
    @test all(rp.JumpProcess.X.p .== [rp.species_to_state[S],
                                      rp.species_to_state[E],
                                      rp.species_to_state[SE],
                                      rp.species_to_state[P]])
end

@testset "Langevin Process - Tests" begin
    # Birth Death Process
    BD = @reaction_network begin
        1.0, ∅ --> A
        0.01, 2A --> A
    end
    A = species(BD)[1]
    A0 = Dict(A => 2.0)
    A_scale = Dict(A => 10.0)
        
    rp = LangevinProcess(BD)
    x = rp.DiffusionProcess.x
    S = prodstoichmat(BD) - substoichmat(BD)
    @test all(rp.DiffusionProcess.f == S*[1.0; 0.01*x[1]*(x[1]-1)/2])
    @test all(rp.DiffusionProcess.σ .== reshape([1.0 + 0.01*x[1]*(x[1]-1)/2],1,1))
    @test all(rp.DiffusionProcess.X.p == [x[1]])

    rp, A0_scaled = langevin_process_setup(BD, A0; scales = A_scale)
    x = rp.DiffusionProcess.x
    @test A0_scaled[1] == A0[A]/A_scale[A]
    @test all(rp.DiffusionProcess.f .≈ S*[1.0/A_scale[A] ; 0.01*x[1]*A_scale[A]*(A_scale[A]*x[1]-1)/2/A_scale[A] ] )
    @test all(rp.DiffusionProcess.σ .≈ reshape([1.0/A_scale[A]^2 + 0.01*x[1]*A_scale[A]*(x[1]*A_scale[A]-1)/2/A_scale[A]^2],1,1))
    @test all(rp.DiffusionProcess.X.p .≈ [A_scale[A]*x[1]])

    # Michaelis-Menten System
    MM = @reaction_network begin
        10.0, S + E --> SE
        0.2, SE --> S + E
        4.0, SE --> P + E
        1.3, P --> S
    end
    S, E, SE, P = species(MM)
    x0 = Dict(S => 10, E => 5, SE => 0, P => 0)
    x_scale = Dict(S => 5, E => 6, SE => 7, P => 8)

    rp, x0_scaled = langevin_process_setup(MM, x0; scales = x_scale)
    x = rp.DiffusionProcess.x
    Smat = (prodstoichmat(MM) - substoichmat(MM))[3:4,:]
    @test all(x0_scaled .== [0.0,0.0])
    ar = [10.0*rp.species_to_state[S]*rp.species_to_state[E];
        0.2*rp.species_to_state[SE];
        4.0*rp.species_to_state[SE];
        1.3*rp.species_to_state[P]]
    @test all(rp.DiffusionProcess.f .≈ Smat*ar ./ [x_scale[SE]; x_scale[P]])
    @test all(rp.DiffusionProcess.σ .≈ [sum(Smat[i,k]*Smat[j,k]*ar[k] for k in 1:4) for i in 1:2, j in 1:2] ./ ([x_scale[SE]; x_scale[P]]*  [x_scale[SE]; x_scale[P]]') )
    @test all(rp.DiffusionProcess.X.p .≈ [rp.species_to_state[S],
                                    rp.species_to_state[E],
                                    rp.species_to_state[SE],
                                    rp.species_to_state[P]])
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