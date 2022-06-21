using Symbolics, LinearAlgebra, MarkovBounds

# initialization via Symbolics.jl
@variables x, y, u, v
a = [0.2*(x-1.2*v)^2*u; 0.2*y + u]
h = [[x + 1, y*v], [0.7*x, (x+1.2*y)^2 + u]]
X = [1.2*(x^2 + 0.2*y)^2 >= 1, (0.2*x - 0.5*y)^2 <= 2.0*x^3]
jp1 = JumpProcess([x,y], a, h, X, controls = [u,v])

@variables x[1:2], u[1:2]
a = [0.2*(x[1]-1.2*u[2])^2*u[1]; 0.2*x[2] + u[1]]
h = [[x[1] + 1, x[2]*u[2]], [0.7*x[1], (x[1]+1.2*x[2])^2 + u[1]]]
X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
jp2 = JumpProcess(x, a, h, X, controls = u)

@variables x[1:2]
a = [0.2*(x[1]-1.2)^2; 0.2*x[2]]
h = [[x[1] + 1, x[2]], [0.7*x[1], (x[1]+1.2*x[2])^2]]
X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
jp3 = JumpProcess(x, a, h, X)

@variables x[1:2], u[1:2]
x = collect(x)
u = collect(u)
a = [0.2*(x[1]-1.2*u[2])^2*u[1]; 0.2*x[2] + u[1]]
h = [[x[1] + 1, x[2]*u[2]], [0.7*x[1], (x[1]+1.2*x[2])^2 + u[1]]]
X = [1.2*(x[1]^2 + 0.2*x[2])^2 >= 1, (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3]
jp4 = JumpProcess(x, a, h, X, controls = u)

# specification via DynamicPolynomials
@polyvar x[1:2]
@polyvar u[1:2]
a = [0.2*(x[1]-1.2*u[2])^2*u[1]; 0.2*x[2] + u[1]]
h = [[x[1] + 1, x[2]*u[2]], [0.7*x[1], (x[1]+1.2*x[2])^2 + u[1]]]
X = @set(1.2*(x[1]^2 + 0.2*x[2])^2 >= 1 && (0.2*x[1] - 0.5*x[2])^2 <= 2.0*x[1]^3)
jp5 = JumpProcess(x, a, h, X, controls = u)
