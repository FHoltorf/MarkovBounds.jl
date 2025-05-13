module MarkovBounds

using Reexport
using SumOfSquares
using ConcreteStructs

# SumOfSquares and DynamicPolynomials are the main modeling tools to set up
# moment bounding problems --> Reexport those
@reexport using SumOfSquares

@reexport using DynamicPolynomials

using Catalyst: species, speciesmap, reactions, paramsmap, prodstoichmat, substoichmat, ReactionSystem, Reaction

using Symbolics: Num, expand, Arr, ArrayShapeCtx, unwrap, istree, PolyForm
using SymbolicUtils: Pow, Mul, Add, Sym, Term, arguments, BasicSymbolic, Symbolic
using Graphs: edges, vertices, SimpleGraph, add_edge!, add_vertex!, Edge
using MetaGraphs: MetaGraph, MetaDiGraph, props, set_prop!
using Base.Iterators: product

import LinearAlgebra: qr, nullspace, diag, Diagonal, tr, transpose
import Base: show
import Parameters: @unpack
import MultivariatePolynomials.polynomial
import SumOfSquares.SemialgebraicSets.inequalities

# Default types (there ought to be a better way to handle this)
const DynPoly = DynamicPolynomials
const COMMUTATIVE = DynPoly.Commutative{DynPoly.CreationOrder}
const POLYORDER = Graded{LexOrder}
const COEFFTYPE = Float64 

const APL = AbstractPolynomialLike
const POLY = DynPoly.Polynomial{COMMUTATIVE, POLYORDER, COEFFTYPE}
const VAR = DynPoly.Variable{COMMUTATIVE, POLYORDER}
const MONO = DynPoly.Monomial{COMMUTATIVE, POLYORDER}

include("distributed.jl")
include("processes.jl")
include("symbolics.jl")
include("print_routines.jl")
include("utils.jl")
include("distributed_utils.jl")
include("problem_setup.jl")
include("sos_programs_stationary.jl")
include("sos_programs_transient.jl")
include("sos_programs_control.jl")
end
