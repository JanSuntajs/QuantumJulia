# ## Time evolution following ED

using DrWatson
@quickactivate :QuantumJulia

using QuantumJulia
using BenchmarkTools
using Random: rand

# We now select some sensible model parameters.
# Feel free to play around with them.
L::Int64= 10
nup::Int64 = 5
J::Float64 = 1.
Δ::Float64 = 1.
W::Float64 = 5.
fields::Vector{Float64} = W * rand(L)

basis = Basis(L, nup)
# This will create the Hamiltonian, but not build it yet:
ham = XXZHamiltonian(basis, J, Δ, fields)

@btime buildham!(ham)

println(ham._isset)

# ## Diagonalization
# We now use the full exact diagonalization to obtain
# all eigenvalues and eigenvectors of the Hamiltonian.
# Proceeed with care as this can be very resource intensive.
λ, V = diagonalize(ham; method=:dense);

# After diagonalization, we can now proceed with the time evolution,
# and, for instance, calculate the survival probability. 
# We will be interested in calculating the following quantity:
# $$
#   P(t) = |\langle \psi(0) | \psi(t) \rangle|^2
# $$
#
# Expanding, we need to evaluate the following:
# $$
#  P(t) = |\sum_\alpha |c_\alpha|^2 \exp(-i E_\alpha t)|^2.
# $$
# Each $c_\alpha$ is just the result of a dot product between the
# column of the V matrix and the initial state vector. We have
# conveniently prepared all the neccessary bits in the SurvivalProbability
# struct. To see the implementation, see `src/states.jl`.
# 
# First, we prepare the random Gaussian state:
ψ0 = random_gaussian_state(basis.nstates)
sprob = SurvivalProbability(λ, V, ψ0)
# We have now just initiated the survival probability struct but
# have not yet calculated the survival probability at a given time.


