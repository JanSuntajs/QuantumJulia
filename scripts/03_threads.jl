# This script demonstrates the use of multithreading in Julia to apply a quantum operator
# to a state vector in a thread-safe manner. It is recommended to run this script as a 
# standalone file rather than in a notebook, as multithreading behavior can sometimes be 
# unclear in notebook environments.
#
# # to run, do
#
# ```bash
# julia --threads=auto --project scripts/03_threads.jl
# ```
#

using DrWatson
@quickactivate :QuantumJulia

using BenchmarkTools
import QuantumJulia.Operators as ops
import LinearAlgebra: norm

"""
    apply_sy!(ψ0::Vector{T}, ψ1::Vector{T}, site::N) where {T <: Complex, N <: Integer}

Apply the `sy` operator to the state vector `ψ0` and update `ψ1` in a thread-safe manner.

# Arguments
- `ψ0::Vector{T}`: Input state vector.
- `ψ1::Vector{T}`: Output state vector to be updated.
- `site::N`: The site on which the `sy` operator acts.

# Notes
To avoid race conditions, `ψ0` is left unchanged and only `ψ1` is updated.
"""
@inline function apply_sy!(ψ0::Vector{T}, 
    ψ1::Vector{T}, site::N) where {T <: Complex, N <: Integer}

    Threads.@threads for i::N in eachindex(ψ1)

        # we have to take i-1, the actual
        # states go from 0 to 2^N-1
        factor, oldstate = ops.sy(i-1, site)
        # this is threadsafe - only a given
        # index i accesses ψ1[i] while
        # ψ0 is not changed; basically, we reverse
        # the approach of the Hamiltonian construction.
        # rather than looking which states can be accessed
        # from the original state, we do the reverse - 
        # we look which states are connected to the final one
        ψ1[i] += factor * ψ0[oldstate + 1]
    end  # i loop

end

# This should give inconsistent results if run in parallel
# due to race conditions
@inline function apply_sy_wrong(ψ0::Vector{T}, 
    ψ1::Vector{T}, site::N) where {T <: Complex, N <: Integer}
    Threads.@threads for i::N in eachindex(ψ1)

        # we have to take i-1, the actual
        # states go from 0 to 2^N-1
        factor, newstate = ops.sy(i-1, site)

        # here, we are appending to the index
        # newstate + 1, which is not threadsafe;
        # many processes may be trying to access
        # it at once.
        ψ1[newstate + 1] += factor * ψ0[i]
    end  # i loop

end
"""
    make_state(N)

Create a normalized quantum state vector of size `2^N`.

# Arguments
- `N`: The number of qubits.

# Returns
- A normalized state vector of type `Vector{ComplexF64}`.
"""
@inline function make_state(N)
    ψ = ones(ComplexF64, 2^N)
    return ψ / norm(ψ)
end

# Example of parallelization
N = 14
ψ0 = make_state(N)
ψ1 = similar(ψ0)



# Benchmark the apply_sy! function
@btime apply_sy!(ψ0, ψ1, 1) setup=(ψ0 = make_state(N); ψ1 = similar(ψ0))

