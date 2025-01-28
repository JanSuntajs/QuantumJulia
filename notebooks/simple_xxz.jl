

using DrWatson
@quickactivate :QuantumJulia

using QuantumJulia
using BenchmarkTools
using KrylovKit

# Feel free to update and ramp these parameters up
L = 12
nup = 6
J = 1.
Δ = 1.

const basis = Basis(L, nup)

ham = XXZHamiltonian(basis, J, Δ)

@btime buildham!(ham)

println(ham._isset)

# Diagonalize
# we use the method eigsolve from KrylovKit in this
# case. We use the default settings of the solver, but the
# formulation of diagonalize allows for both the specification of
# the additional positional and keyword arguments:
# [KrylovKit](https://jutho.github.io/KrylovKit.jl/stable/man/eig/)
λ, V = diagonalize(ham; method=:krylov)

# find the ground state
println(λ)
