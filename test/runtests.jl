using Powersense
using GLPK
using Ipopt
using Test



# @testset "Unit tests" begin
#     include("unittests.jl")
# end


@testset "MathOptInterface" begin
    include("opt/MOI_wrapper.jl")
end

@testset "External Solver Attributes Implementation with Toy Example" begin
    include("opt/ext_solver.jl")
    @test typeof(optimizer_solver) == MOI.OptimizerWithAttributes
    @test isapprox(xsol, -1.0, rtol=1e-4)
    @test isapprox(ysol, -1.0, rtol=1e-4)
    @test status == MOI.LOCALLY_SOLVED
    include("../examples/opt/opt_example.jl")
    @test isapprox(xsol, -1.0, rtol=1e-4)
    @test isapprox(ysol, -1.0, rtol=1e-4)
end

@testset "ACOPF Formulations" begin
    include("../examples/opf/ACOPF_formulations_example.jl")
    include("opf/ACOPF_formulations.jl")
end

