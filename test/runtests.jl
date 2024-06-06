using QEDcore
using Test

@testset "QEDcore.jl" begin
    @testset "dummy" begin
        @test QEDcore.greet("World!") == "Hello World!"
    end
end
