using QEDcore
using Test
using SafeTestsets

begin
    @time @safetestset "four momentum" begin
        include("four_momentum.jl")
    end

    @time @safetestset "gamma matrices" begin
        include("gamma_matrices.jl")
    end

    @time @safetestset "Lorentz vector" begin
        include("lorentz_vector.jl")
    end
end
