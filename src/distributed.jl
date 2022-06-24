export Partition, Singleton, @singleton

mutable struct Singleton <: AbstractSemialgebraicSet
    x::AbstractVector
end
Singleton(x::Real) = Singleton([x])
Singleton(x::Tuple) = Singleton([x...])
macro singleton(x)
    return :(Singleton($(esc(x))))
end

struct Partition
    graph
    get_vertex
end