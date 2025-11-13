using Test
using SafeTestsets

begin
    @safetestset "Ideal system" begin
        include("electron_system/ideal.jl")
    end
    @safetestset "Screening" begin
        include("electron_system/screening.jl")
    end
    @safetestset "Interacting system" begin
        include("electron_system/interacting.jl")
    end
end
