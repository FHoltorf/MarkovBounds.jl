Base.show(io::IO, MP::JumpProcess) = println(io,
"Jump Process
-----------------------
  States: $(print_wo_type(MP.x))
  Propensities: $(print_wo_type(MP.a))
  Jumps: ["*prod(print_wo_type(MP.h[i])*(i < length(MP.h) ? ", " : "" ) for i in 1:length(MP.h))*"]
  State Space: $(MP.X)")

Base.show(io::IO, RP::ReactionProcess) = println(io,
"Reaction Process is analogous to following $(RP.JumpProcess)")

Base.show(io::IO, MP::DiffusionProcess) = println(io,
"Diffusion Process
-----------------------
  States: $(print_wo_type(MP.x))
  Drift: $(print_wo_type(MP.f))
  Diffusion: $(print_wo_type(MP.σ))
  State space: $(MP.X)")

Base.show(io::IO, MP::JumpDiffusionProcess) = println(io,
"Jump-Diffusion Process
-----------------------
  States: $(print_wo_type(MP.x))
  Propensities: $(print_wo_type(MP.a))
  Jumps: ["*prod(print_wo_type(MP.h[i])*(i < length(MP.h) ? ", " : "" ) for i in 1:length(MP.h))*"]
  Drift: $(print_wo_type(MP.f))
  Diffusion: $(print_wo_type(MP.σ))
  State Space: $(MP.X)")

Base.show(io::IO, CP::ControlProcess) = println(io,
"Control Process
-----------------------
$(CP.MP)-----------------------
Admissible Controls: $(CP.U)
Control Horizon: [0, $(CP.T)]
Objective Function: $(CP.Objective)")

Base.show(io::IO, obj::LagrangeMayer) = println(io,
"𝔼[∫ $(obj.l) dt + $(obj.m)]")

Base.show(io::IO, S::Singleton) = println(io,"{ $(S.x) }")

function print_wo_type(A::AbstractArray)
  str = string(A)
  idx = findfirst('[', string(A))
  return str[idx:end]
end
