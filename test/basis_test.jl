"""
In this test module, we test the correctness
of the implementation of the basis struct. In testing,
we should strive for the tests to be as exhaustive as
possible. The test names and test cases should also be rather
descriptive and self-explanatory so that one knows where a given
test has failed.


"""
module BasisTest

using Test
using QuantumJulia

# Test the Basis struct
@testset "Test Basis implementation" begin

    # Test for full Hilbert space
    # we test whether correct states are
    # generated for a manageable system size
    size = 2
    basis = Basis(size, nothing)
    @test basis.states == [0, 1, 2, 3]

    # in the next tests, we check whether also
    # the spin conservation sectors work as expected.
    size = 2
    basis = Basis(size, 1)
    @test basis.states == [1, 2]

    basis = Basis(size, [1, 2])
    @test basis.states == [1, 2, 3]
    @test basis.nstates == 3

    # slightly larger system, we check whether
    # the number of states is correct
    size = 4
    basis = Basis(size, nothing)
    @test length(basis.states) == 2^size
    @test basis.nstates == 2^size

    # Test for particle conserving sectors
    nup = 2
    basis = Basis(size, nup)
    @test length(basis.states) == binomial(size, nup)
    @test basis.nstates == binomial(size, nup)

    # Test for multiple particle conserving sectors
    nup = [1, 2]
    basis = Basis(size, nup)
    expected_states = sum(binomial(size, nu) for nu in nup)
    @test length(basis.states) == expected_states
    @test basis.nstates == expected_states
end

end # module BasisTest
