"""
In this script, we benchmark the performance of the `buildham!` function.
To run as a script, do the following:

```bash
julia --project=. scripts/01_simple_xxz.jl
```
Here, we have assumed you are running the code from the root folder of the
project. Since we are managing the project using DrWatson, this is not strictly
needed. Even running from elsewhere (and without the --project=.) line will be fine
as the `@quickactivate` macro will take care of the environment activation.

"""

using DrWatson
@quickactivate :QuantumJulia

using QuantumJulia
using BenchmarkTools
using Random: rand

# Feel free to update and ramp these parameters up
const L::Int64 = 16
const nup::Int64 = 8
const J::Float64 = 1.
const Δ::Float64 = 1.
const W::Float64 = 5.
const fields::Vector{Float64} = W * rand(L)

const basis::Basis = Basis(L, nup)
# This will create the Hamiltonian, but not build it yet.
ham::XXZHamiltonian = XXZHamiltonian(basis, J, Δ, fields)

@btime buildham!(ham)

println(ham._isset)

# Diagonalize
# we use the method eigsolve from KrylovKit in this
# case. We use the default settings of the solver, but the
# formulation of diagonalize allows for both the specification of
# the additional positional and keyword arguments
# https://jutho.github.io/KrylovKit.jl/stable/man/eig/
println("Starting diagonalization!")
@btime λ, V = diagonalize(ham; method=:dense);
println("Diagonalization done!")
λ, V = diagonalize(ham; method=:dense);
# Now calculate the survival probability.
ψ0 = random_gaussian_state(basis.nstates)
@time sprob = SurvivalProbability(λ, V, ψ0)

times = exp10.(range(-3, stop=4, length=1000));
probvals = similar(times)
@time probvals .= sprob.(times);
