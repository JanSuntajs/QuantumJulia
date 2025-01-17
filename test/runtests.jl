using DrWatson
@quickactivate :QuantumJulia

using Test

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
println("Starting tests")

@testset "QuantumJulia tests" begin

    @testset "Basis tests" begin
        include("./basis_test.jl")
    end

end
