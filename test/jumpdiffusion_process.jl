using Symbolics, LinearAlgebra, MarkovBounds

# initialization via Symbolics.jl
@variables x, y, u, v
f = [1.5*x*(x-1)^2*u-v^2; 1.2*y+(u+0.2)^2]
σ = Matrix(Diagonal([0.2*x^2; 0.2*y]))
a = [0.2*(x-1.2*v)^2*u; 0.2*y + u]
h = [[x + 1, y*v], [0.7*x, (x+1.2*y)^2 + u]]
X = [1.2*(x^2 + 0.2*y)^2 >= 1, (0.2*x - 0.5*y)^2 <= 2.0*x^3]
jpd1 = JumpDiffusionProcess([x,y], a, h, f, σ, X, controls = [u,v])

@variables x[1:2], u[1:2]
f = [1.5*x[1]*(x[1]-1)^2*u[1]-u[2]^2; 1.2*x[2]+(u[1]+0.2)^2]
σ = Matrix(Diagonal([0.2*x[1]^2; 0.2*x[2]]))
a = [0.2*(x[1]-1.2*u[2])^2*u[1]; 0.2*x[2] + u[1]]
h = [[x[1] + 1, x[2]*u[2]], [0.7*x[1], (x[1]+1.2*x[2])^2 + u[1]]]
X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
jpd2 = JumpDiffusionProcess(x, a, h, f, σ, X, controls = u)

@variables x[1:2]
f = [1.5*x[1]*(x[1]-1)^2; 1.2*x[2]]
σ = Matrix(Diagonal([0.2*x[1]^2; 0.2*x[2]]))
a = [0.2*(x[1]-1.2)^2; 0.2*x[2]]
h = [[x[1] + 1, x[2]], [0.7*x[1], (x[1]+1.2*x[2])^2]]
X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
jpd3 = JumpDiffusionProcess(x, a, h, f, σ, X)

@variables x[1:2], u[1:2]
x = collect(x)
u = collect(u)
f = [1.5*x[1]*(x[1]-1)^2*u[1]-u[2]^2; 1.2*x[2]+(u[1]+0.2)^2]
σ = Matrix(Diagonal([0.2*x[1]^2; 0.2*x[2]]))
a = [0.2*(x[1]-1.2*u[2])^2*u[1]; 0.27*x[2] + u[1]]
h = [[x[1] + 1, x[2]*u[2]], [0.7*x[1], (x[1]+1.2*x[2])^2 + u[1]]]
X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
jpd4 = JumpDiffusionProcess(x, a, h, f, σ, X, controls = u)

# specification via DynamicPolynomials
@polyvar x[1:2]
@polyvar u[1:2]
f = [1.5*x[1]*(x[1]-1)^2*u[1]-u[2]^2; 1.2*x[2]+(u[1]+0.2)^2]
σ = Matrix(Diagonal([0.2*x[1]^2; 0.2*x[2]]))
a = [0.2*(x[1]-1.2*u[2])^2*u[1]; 0.2*x[2] + u[1]]
h = [[x[1] + 1, x[2]*u[2]], [0.7*x[1], (x[1]+1.2*x[2])^2 + u[1]]]
X = @set(1.2*(x[1]^2 + 0.2*x[2])^2 >= 1 && (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3)
jpd5 = JumpDiffusionProcess(x, a, h, f, σ, X, controls = u)
