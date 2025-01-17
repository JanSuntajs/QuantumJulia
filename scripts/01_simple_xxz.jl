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
using KrylovKit

# feel free to update and ramp these up
const L = 12
const nup = 6
const J = 1.
const Δ = 1.

const basis = Basis(L, nup)

ham = XXZHamiltonian(basis, J, Δ)

@btime buildham!(ham)

println(ham._isset)

# diagonalize
# we use the method eigsolve from KrylovKit in this
# case. We use the default settings of the solver, but the
# formulation of diagonalize allows for both the specification of
# the additional positional and keyword arguments
# https://jutho.github.io/KrylovKit.jl/stable/man/eig/
λ, V = diagonalize(ham; method=:krylov)

println(λ)

# find the ground state