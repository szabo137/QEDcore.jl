using QEDcore
using Test
using SafeTestsets

begin
    # interfaces
    @time @safetestset "process interface" begin
        include("interfaces/process.jl")
    end

    @time @safetestset "computation setup interface" begin
        include("interfaces/setup.jl")
    end

    @time @safetestset "four momentum" begin
        include("four_momentum.jl")
    end

    @time @safetestset "gamma matrices" begin
        include("gamma_matrices.jl")
    end

    @time @safetestset "Lorentz vector" begin
        include("lorentz_vector.jl")
    end

    @time @safetestset "Dirac tensors" begin
        include("dirac_tensor.jl")
    end

    @time @safetestset "particle spinors" begin
        include("particle_spinors.jl")
    end

    @time @safetestset "particle spinors" begin
        include("particle_base_states.jl")
    end
end
