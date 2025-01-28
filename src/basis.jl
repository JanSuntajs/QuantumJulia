"""
A struct representing a basis for a quantum system.

# Fields
- `size::T`: Size of the system.
- `nup::Union{T, Vector{T}, Nothing}`: Number of up spins or a vector of up spins.
- `states::Vector{T}`: Vector of basis states.
- `nstates::T`: Number of basis states.

# Constructor
- `Basis(size::T, nup::Union{T, Vector{T}, Nothing})`: Creates a new `Basis` instance.
"""
struct Basis{T<:Integer} 
    size::T
    nup::Union{T, Vector{T}, Nothing}
    states::Vector{T}
    nstates::T

    function Basis(size::T, nup::Union{T, Vector{T}, Nothing}) where {T<:Integer}
        states = _getstates(size, nup)
        nstates = length(states)
        new{T}(size, nup, states, nstates)
    end
end


"""Set state indices for full Hilbert space."""
function _getstates(size::T, nup::Nothing)::Vector{T} where {
        T <: Integer}

    states = collect(T, 0:1<<size-1)


    return states

end

"""Set state indices for particle conserving sectors"""
function _getstates(size::T,
    nup::Union{T, Vector{T}})::Vector{T} where {
        T<:Integer}

    nstates = T(sum(binomial(size, nu) for nu in nup))
    states = Vector{T}(undef, nstates)
    start::T = 0

    for nu in nup
        states_ = selectstates(size, nu)
        states[start+1:start+length(states_)] .= states_
        start += T(binomial(size, nu))
    end

    return states
end