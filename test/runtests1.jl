push!(LOAD_PATH, "/home/projects/projects/nonlinear_optimization/ActiveSetMethodss/src");
using ActiveSetMethods
using GLPK
using Test



# @testset "Unit tests" begin
#     include("unittests.jl")
# end


@testset "MathOptInterface" begin
    include("MOI_wrapper.jl")
end

#=

@testset "External Solver Attributes Implementation with Toy Example" begin
    include("ext_solver.jl")
    @test typeof(optimizer_solver) == MOI.OptimizerWithAttributes
    @test isapprox(xsol, -1.0, rtol=1e-4)
    @test isapprox(ysol, -1.0, rtol=1e-4)
    @test status == MOI.LOCALLY_SOLVED
end

@testset "toy_example.jl" begin
    include("../examples/toy_example.jl")
    @test isapprox(xsol, -1.0, rtol=1e-4)
    @test isapprox(ysol, -1.0, rtol=1e-4)
end

@testset "opf.jl" begin
    include("../examples/acopf/opf.jl")
    run_opf("../examples/acopf/case3.m")
end 
=#
