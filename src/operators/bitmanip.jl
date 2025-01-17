"""
This module contains routines for bit manipulations which
are applicable to both hard-core bosonic (spin 1/2) or
fermionic systems -> to all systems with two local degrees
of freedom.

Routines defined here are:

countSetbits(n) - returns the number of set (nonzero)
                  bits in a state. Even though one
                  could technically use bin(n).count('1')
                  to achieve the same result, this approach
                  cannot be used in the nopython mode of
                  numba jit.

select_states(L, nup) - a routine for selecting the states
                        with an appropriate number of particles
                        or with an appropriate projection of
                        the total spin along the z-axis.

bitflip(state, bit) - a routine for flipping a bit in a state
                      at a position 'bit' into the opposite
                      state.

getbit(state, bit) - get the value of a bit in a state at a
                     position 'bit'.
"""


"""
    bitflip(state::T, bit::N) -> T

Flip a bit at position `bit` of a state `state`.

# Arguments
- `state::T where T <: Unsigned`: state on which we act.
- `bit::N where N <: Unsigned`: bit which we flip.

# Returns
- `T`: state with a flipped `bit`.

# Example
```julia
julia> bitflip(0b0101, 1)
0b0100
```
"""
@inline function bitflip(state::T, bit::N)::T where {T <: Integer, N <: Integer}
    # Use XOR to flip the bit at the specified position
    return state ⊻ (1 << bit)
end


"""
    getbit(state::T, bit::N) -> N

Get the value of a `bit` in a state `state`.

# Arguments
- `state::T where T <: Unsigned`: state on which we act.
- `bit::N where N <: Unsigned`: the checked bit.

# Returns
- `N`: the `bit` value of a `state`.

# Example
```julia
julia> getbit(0b0101, 2)
1
```
"""
@inline function getbit(state::T, bit::N)::N where {T <: Integer, N <: Integer}
    # Shift the state right by `bit` positions and mask with 1 to get the bit value
    return ((state >> bit) & 1)
end


"""
    countsetbits(state::T) -> T

Count the number of set bits (ones) in a given `state`.

# Arguments
- `state::T where T <: Unsigned`: the state we are checking.

# Returns
- `T`: the number of set bits.

# Example
```julia
julia> countsetbits(0b0101)
2
```
"""
@inline function countsetbits(state::T)::T where {T<:Integer}
    # Use the built-in function to count the number of set bits
    return count_ones(state)
end



"""
    selectstates(L::T, nup::T) -> Vector{T}

Select (basis) states with a given number of particles.
This will be applicable in calculations with large matrices
and on-the-fly calculations.

# Arguments
- `L::T where T <: Unsigned`: system size.
- `nup::T where T <: Unsigned`: number of up spins.

# Returns
- `Vector{T}`: vector of selected states.

# Example
```julia
julia> selectstates(4, 2)
4-element Vector{UInt64}:
 0x0000000000000003
 0x0000000000000005
 0x0000000000000006
 0x0000000000000009
```
"""
@inline function selectstates(L::T, nup::T)::Vector{T} where {T<:Integer}
    @assert nup <= L
    # Preallocate the vector for selected states
    states = Vector{T}(undef, binomial(L, nup))

    idx_in::T = 0
    @inbounds for i::T in 0:1<<L - 1
        # Check if the number of set bits equals nup
        if countsetbits(i) == nup
            idx_in += 1
            states[idx_in] = i
        end
    end
    return states
end