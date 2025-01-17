"""
Implementation of spin 1/2
operators to help in Hamiltonian
creation procedure.

"""

"""
    sx(state::T, bit::N) -> Tuple{Float64,T}

Spin-1/2 Pauli-X operator. Flips the bit at the specified position.

# Arguments
- `state::T`: The quantum state.
- `bit::N`: The bit position to flip.

# Returns
- `Tuple{Float64,T}`: A tuple containing the factor (0.5) and the new state with the bit flipped.

# Example
```julia
julia> sx(0b0101, 1)
(0.5, 0b0100)
```
"""
@inline function sx(state::T, bit::N)::Tuple{Float64,T} where {T <: Integer, N<:Integer}

    return 0.5, bitflip(state, bit)
end

"""
    sy(state::T, bit::N) -> Tuple{ComplexF64,T}

Spin-1/2 Pauli-Y operator. Flips the bit at the specified position and applies a phase factor.

# Arguments
- `state::T`: The quantum state.
- `bit::N`: The bit position to flip.

# Returns
- `Tuple{ComplexF64,T}`: A tuple containing the phase factor and the new state with the bit flipped.

# Example
```julia
julia> sy(0b0101, 1)
(0.0 + 0.5im, 0b0100)
```
"""
@inline function sy(state::T, bit::N)::Tuple{ComplexF64,T} where {T <: Integer, N<:Integer}

    factor = 0.5im * (-1)^(1 - getbit(state, bit))

    return factor, bitflip(state, bit)
end

"""
    sz(state::T, bit::N) -> Tuple{Float64,T}

Spin-1/2 Pauli-Z operator. Applies a phase factor based on the bit at the specified position.

# Arguments
- `state::T`: The quantum state.
- `bit::N`: The bit position to check.

# Returns
- `Tuple{Float64,T}`: A tuple containing the phase factor and the unchanged state.

# Example
```julia
julia> sz(0b0101, 2)
(-0.5, 0b0101)
```
"""
@inline function sz(state::T, bit::N)::Tuple{Float64,T} where {T <: Integer, N<:Integer}

    factor = 0.5 * (-1.0)^(1 - getbit(state, bit))

    return factor, state
end

"""
    id2(state::T, bit::N) -> Tuple{Float64,T}

Identity operator for spin-1/2. Returns the state unchanged.

# Arguments
- `state::T`: The quantum state.
- `bit::N`: The bit position (unused).

# Returns
- `Tuple{Float64,T}`: A tuple containing the factor (1.0) and the unchanged state.

# Example
```julia
julia> id2(0b0101, 2)
(1.0, 0b0101)
```
"""
@inline function id2(state::T, bit::N)::Tuple{Float64,T} where {T <: Integer, N<:Integer}

    return 1.0, state
end

"""
    sp(state::T, bit::N) -> Tuple{Float64,T}

Spin-1/2 raising operator. Raises the bit at the specified position.

# Arguments
- `state::T`: The quantum state.
- `bit::N`: The bit position to raise.

# Returns
- `Tuple{Float64,T}`: A tuple containing the factor and the new state with the bit raised.

# Example
```julia
julia> sp(0b0101, 1)
(1.0, 0b0100)
```
"""
@inline function sp(state::T, bit::N)::Tuple{Float64,T} where {T <: Integer, N<:Integer}

    bitval = getbit(state, bit)

    if !iszero(bitval)
        return 1 - bitval, state
    else
        return 1 - bitval, bitflip(state, bit)
    end
end

"""
    sm(state::T, bit::N) -> Tuple{Float64,T}

Spin-1/2 lowering operator. Lowers the bit at the specified position.

# Arguments
- `state::T`: The quantum state.
- `bit::N`: The bit position to lower.

# Returns
- `Tuple{Float64,T}`: A tuple containing the factor and the new state with the bit lowered.

# Example
```julia
julia> sm(0b0101, 1)
(0.0, 0b0101)
```
"""
@inline function sm(state::T, bit::N)::Tuple{Float64,T} where {T <: Integer, N<:Integer}

    bitval = getbit(state, bit)

    if !iszero(bitval)
        return bitval, bitflip(state, bit)
    else
        return bitval, state
    end
end