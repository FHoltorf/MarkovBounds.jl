export Partition, Singleton, @singleton, SemidiscreteSet

mutable struct Singleton
    x::AbstractVector
end
Singleton(x::Real) = Singleton([x])
Singleton(x::Tuple) = Singleton([x...])
macro singleton(x)
    return :(Singleton($(esc(x))))
end

mutable struct SemidiscreteSet
    discrete # vector of pairs x => [vals] (each referencing a singleton in the projected space)
    continuous # BasicSemialgebraicSet
end

struct Partition
    graph
    get_vertex
end