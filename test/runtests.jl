using QEDcore
using Test
using SafeTestsets

begin
    # algebraic objects
    @time @safetestset "four momentum" begin
        include("algebraic_objects/four_momentum.jl")
    end

    @time @safetestset "gamma matrices" begin
        include("algebraic_objects/gamma_matrices.jl")
    end

    @time @safetestset "Lorentz vector" begin
        include("algebraic_objects/lorentz_vector.jl")
    end

    @time @safetestset "Dirac tensors" begin
        include("algebraic_objects/dirac_tensor.jl")
    end

    # particles
    @time @safetestset "particle spinors" begin
        include("particles/spinors.jl")
    end
end
