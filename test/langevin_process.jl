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