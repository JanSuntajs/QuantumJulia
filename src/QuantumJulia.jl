"""
    module QuantumJulia

This module serves as a pedagogical introduction to Julia for quantum 
Hamiltonian construction and manipulation. It includes the necessary 
functions and structures to define and build Hamiltonians for spin-1/2 
systems.

# Explanation
In Julia, a module is a way to group related functions, types, and 
variables together. It helps in organizing code and managing namespaces. 
Modules can be nested, and each module has its own global scope.

# DrWatson
DrWatson is a Julia package that provides tools for scientific project 
management. It helps in organizing code, data, and results in a 
reproducible manner. The `@quickactivate` macro is used to quickly 
activate the project environment.

# Disclaimer
This module and its contents were generated using GitHub Copilot, an AI 
programming assistant. The code and explanations provided are intended 
for educational purposes.

# Usage
To use this module, simply include it in your Julia script or REPL 
session:
```julia
using QuantumJulia
```
"""
module QuantumJulia

    using DrWatson
    @quickactivate :QuantumJulia

    using Distributions
    using Random
    using LinearAlgebra
    using LinearAlgebra: eigen, eigen!
    using KrylovKit: eigsolve
    using SparseArrays

    include("./operators/operators.jl")
    using .Operators
    include("./basis.jl")
    include("./hamiltonian.jl")
    include("./buildham.jl")
    include("./diagonalize.jl")

    export Basis, XXZHamiltonian, buildham!, diagonalize
end