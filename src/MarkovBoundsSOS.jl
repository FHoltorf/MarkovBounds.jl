module MarkovBoundsSOS

using Reexport
using SumOfSquares

# SumOfSquares and DynamicPolynomials are the main modeling tools to set up
# moment bounding problems --> Reexport those
@reexport using SumOfSquares: @set, SemialgebraicSets

@reexport using DynamicPolynomials: @polyvar, PolyVar, AbstractPolynomialLike,
                                    Polynomial, Term, subs, polynomial, differentiate,
                                    MonomialVector, maxdegree

# MTK/Symbolics/SymbolicUtils/Catalyst can also be used to specify moment problem
# data but not as its main approach
using Catalyst: @reaction_network, @parameters, species, speciesmap,
                reactions, paramsmap, prodstoichmat, substoichmat,
                ReactionSystem, Reaction

using Symbolics: Num, expand
using SymbolicUtils: Add, Mul, Pow, Sym, arguments


import LinearAlgebra: qr, nullspace, diag, Diagonal, tr, transpose
import Base: show
import Parameters: @unpack

const PV = PolyVar
const APL = AbstractPolynomialLike

include("processes.jl")
include("symbolics.jl")
include("SOSPrograms.jl")
include("distributed.jl")
include("print_wrapper.jl")
include("utils.jl")

export JumpProcess, ReactionProcess, DiffusionProcess, JumpDiffusionProcess, ControlProcess,
       LagrangeMayer, Lagrange, Mayer, ExitProbability, TerminalSetProbability,
       setup_reaction_process, inf_generator, extended_inf_generator,
       stationary_polynomial, stationary_mean, stationary_variance, stationary_covariance_ellipsoid,
       transient_polynomial, transient_mean, transient_variance, transient_covariance_ellipsoid,
       optimal_control, value_function
end
