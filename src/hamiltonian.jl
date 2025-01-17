"""
    struct XXZHamiltonian{T <: Number, N <: Integer, M <: Basis{N}}

A mutable struct representing the XXZ Hamiltonian. In this case, it is important that
the struct is mutable, since we will be modifying the `mat` field in-place when creating
the Hamiltonian matrix or modifying other parameters. In cases where we do not change
the struct properties, it is better to use the default (immutable) struct if performance
is critical.

# Fields
- `basis::M`: The basis for the Hamiltonian.
- `mat::SparseMatrixCSC{T, N}`: The Hamiltonian matrix.
- `Delta::T`: The interaction strength.
- `J::T`: The coupling constant.
- `fields::Vector{T}`: The magnetic field values.
- `_full_hilbert_space::Bool`: A flag indicating if the full Hilbert space is used.

# Notes
For those coming from object-oriented languages, note that the concept of a class and related
object-oriented programming concepts are not really present in Julia. Instead, one has structs.
In this particular implementation, we also include a constructor function; one of its functionalities
is to populate the fields with zeros if the user does not provide the `fields` argument.
# Constructor
- `XXZHamiltonian(basis::M, Delta::T, J::T, fields::Vector{T}=zeros(T, basis.size))`: Creates a new XXZHamiltonian instance.
"""
mutable struct XXZHamiltonian{T <: Number, N <: Integer, M <: Basis{N}}

    basis::M
    mat::SparseMatrixCSC{T, N}
    Delta::T
    J::T
    fields::Vector{T}
    _full_hilbert_space::Bool
    _isset::Bool


    function XXZHamiltonian(basis::M, Delta::T, J::T, 
        fields::Vector{T}=zeros(T, basis.size)) where {T <: Number, N <: Integer, M <: Basis{N}}

        @assert length(fields) == basis.size "Length of fields must match the number of sites."
        mat = spzeros(T, N, basis.nstates, basis.nstates)
        _full_hilbert_space = basis.nstates == 2^basis.size
        new{T, N, M}(basis, mat, Delta, J, fields, _full_hilbert_space, false)
    end
end  # struct XXZHamiltonian






