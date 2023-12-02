export Partition, Singleton, @singleton, SemidiscreteSet

@concrete mutable struct Singleton
    x<:AbstractVector
end
Singleton(x::Real) = Singleton([x])
Singleton(x::Tuple) = Singleton(collect(x))
macro singleton(x)
    return :(Singleton($(esc(x))))
end

@concrete mutable struct SemidiscreteSet
    discrete # vector of pairs x => [vals] (each referencing a singleton in the projected space)
    continuous # BasicSemialgebraicSet
end

@concrete struct Partition
    graph
    get_vertex
end