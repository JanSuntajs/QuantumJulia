const _allowed_diagonalization_methods = [:dense, :dense_inplace, :krylov,]



function diagonalize(ham::XXZHamiltonian, args...; method::Symbol=:dense, kwargs...)

    # do some checking
    if method ∉ _allowed_diagonalization_methods
        throw(ArgumentError("Method $method not allowed. Allowed 
        methods are $_allowed_diagonalization_methods"))
    end

    if !ham._isset
        throw(ArgumentError("Hamiltonian matrix not set. Please set 
        the Hamiltonian matrix before diagonalizing."))
    end

    if method == :dense
        mat = issparse(ham.mat) ? Matrix(ham.mat) : ham.mat
        return eigen(mat, args...; kwargs...)
    elseif method == :dense_inplace
        mat = issparse(ham.mat) ? Matrix(ham.mat) : ham.mat
        return eigen!(mat; kwargs...)
    elseif method == :krylov
        # krylov method
        mat = issparse(ham.mat) ? ham.mat : sparse(ham.mat)
        return eigsolve(mat, args...; kwargs...)
    end

end