using MarkovBounds
using Catalyst
using Symbolics
using Test
using SCS # Hypatia 
using LinearAlgebra

solver = optimizer_with_attributes(SCS.Optimizer, "verbose" => false)

include("reaction_process.jl")
include("langevin_process.jl")
include("control_process.jl")
include("jump_process.jl")
include("drift_process.jl")
include("diffusion_process.jl")
include("jumpdiffusion_process.jl")
include("stationary_bounds.jl")
include("control_bounds.jl")