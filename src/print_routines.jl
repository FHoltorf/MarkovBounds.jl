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

Base.show(io::IO, LP::LangevinProcess) = println(io,
"Chemical Langevin Equation is analogous to following $(LP.DiffusionProcess)")

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

Base.show(io::IO, bound::Bound) = println(io,
"$(bound.value) $(interpret_status(termination_status(bound.model), primal_status(bound.model)))")

Base.show(io::IO, S::Singleton) = println(io,"{ $(S.x) }")

function print_wo_type(A::AbstractArray)
  str = string(A)
  idx = findfirst('[', string(A))
  return str[idx:end]
end

function interpret_status(opt_status, feas_status)
  if feas_status == MOI.FEASIBLE_POINT
      status = "valid"
  elseif feas_status == MOI.NEARLY_FEASIBLE_POINT
      status = "likely valid"
  elseif feas_status == MOI.INFEASIBLE_POINT
      status = "invalid"
  else
      status = "unknown"
  end

  if opt_status == MOI.OPTIMAL
      status *= " (optimal)"
  elseif opt_status == MOI.ALMOST_OPTIMAL
      status *= " (near optimal)"
  elseif opt_status == MOI.INFEASIBLE
      status *= " (infeasible)"
  elseif opt_status == MOI.SLOW_PROGRESS
      status *= " (potentially suboptimal)"
  elseif opt_status == MOI.MEMORY_LIMIT
      status *= " (memory limit exceeded)"
  elseif opt_status == MOI.OBJECTIVE_LIMIT
      status *= " (objective limit exceeded)"
  else
      status *= " (unkown)"
  end
  return status
end