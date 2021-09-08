using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

const optimizer = ActiveSetMethods.Optimizer()
MOI.set(optimizer, MOI.RawParameter("external_optimizer"), GLPK.Optimizer)
MOI.set(optimizer, MOI.RawParameter("max_iter"), 1000)

const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4,
                               optimal_status=MOI.LOCALLY_SOLVED)
# DualObjectiveValue is not implemented, so ActiveSetMethods does not pass the tests that
# query it.
# TODO: Consider implementing DualObjectiveValue for purely linear problems.
const config_no_duals = MOIT.TestConfig(atol=1e-4, rtol=1e-4, duals=false,
                                        optimal_status=MOI.LOCALLY_SOLVED)

@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "ActiveSetMethods"
end

@testset "supports_default_copy_to" begin
    @test MOIU.supports_default_copy_to(optimizer, false)
    @test !MOIU.supports_default_copy_to(optimizer, true)
end

@testset "Unit" begin
    bridged = MOIB.full_bridge_optimizer(
        ActiveSetMethods.Optimizer(external_optimizer = GLPK.Optimizer),
        Float64)
    # A number of test cases are excluded because loadfromstring! works only
    # if the solver supports variable and constraint names.
    exclude = ["delete_variable", # Deleting not supported.
               "delete_variables", # Deleting not supported.
               "getvariable", # Variable names not supported.
               "solve_zero_one_with_bounds_1", # Variable names not supported.
               "solve_zero_one_with_bounds_2", # Variable names not supported.
               "solve_zero_one_with_bounds_3", # Variable names not supported.
               "getconstraint", # Constraint names not suported.
               "variablenames", # Variable names not supported.
               "solve_with_upperbound", # loadfromstring!
               "solve_with_lowerbound", # loadfromstring!
               "solve_integer_edge_cases", # loadfromstring!
               "solve_affine_lessthan", # loadfromstring!
               "solve_affine_greaterthan", # loadfromstring!
               "solve_affine_equalto", # loadfromstring!
               "solve_affine_interval", # loadfromstring!
               "get_objective_function", # Function getters not supported.
               "solve_constant_obj",  # loadfromstring!
               "solve_blank_obj", # loadfromstring!
               "solve_singlevariable_obj", # loadfromstring!
               "solve_objbound_edge_cases", # ObjectiveBound not supported.
               "solve_affine_deletion_edge_cases", # Deleting not supported.
               "solve_unbounded_model", # `NORM_LIMIT`
               "number_threads", # NumberOfThreads not supported
               "delete_nonnegative_variables", # get ConstraintFunction n/a.
               "update_dimension_nonnegative_variables", # get ConstraintFunction n/a.
               "delete_soc_variables", # VectorOfVar. in SOC not supported
               "solve_result_index", # DualObjectiveValue not supported
               ]
    MOIT.unittest(bridged, config_no_duals, exclude)
end

@testset "MOI Linear tests" begin
    exclude = ["linear8a", # Behavior in infeasible case doesn't match test.
               "linear12", # Same as above.
               "linear8b", # Behavior in unbounded case doesn't match test.
               "linear8c", # Same as above.
               "linear7",  # VectorAffineFunction not supported.
               "linear15", # VectorAffineFunction not supported.
               ]
    model_for_ActiveSetMethods = MOIU.UniversalFallback(MOIU.Model{Float64}())
    linear_optimizer = MOI.Bridges.Constraint.SplitInterval{Float64}(
                         MOIU.CachingOptimizer(model_for_ActiveSetMethods, optimizer))
    MOIT.contlineartest(linear_optimizer, config_no_duals, exclude)
    # Tests setting bounds of `SingleVariable` constraint
    MOIT.linear4test(optimizer, config_no_duals)
end

MOI.empty!(optimizer)

# This tests qp1, qp2, and qp3.
@testset "MOI QP tests" begin
    qp_optimizer = MOIU.CachingOptimizer(MOIU.Model{Float64}(), optimizer)
    # SLP returns the max iteration limit for these instances.
    exclude = ["qp1", "qp2"]
    MOIT.qptest(qp_optimizer, config_no_duals, exclude)
    # MOIT.qptest(qp_optimizer, config_no_duals)
end

MOI.empty!(optimizer)

@testset "MOI QCQP tests" begin
    qp_optimizer = MOIU.CachingOptimizer(MOIU.Model{Float64}(), optimizer)
    exclude = ["qcp1"] # VectorAffineFunction not supported.
    MOIT.qcptest(qp_optimizer, config_no_duals, exclude)
end

MOI.empty!(optimizer)

@testset "MOI NLP tests" begin
    # TODO: These instances do not converge.
    exclude = ["hs071_no_hessian", "hs071"]
    MOIT.nlptest(optimizer, config, exclude)
    # MOIT.nlptest(optimizer, config)
end

@testset "Testing getters" begin
    MOI.Test.copytest(MOI.instantiate(ActiveSetMethods.Optimizer, with_bridge_type=Float64), MOIU.Model{Float64}())
end

@testset "Bounds set twice" begin
    MOI.Test.set_lower_bound_twice(optimizer, Float64)
    MOI.Test.set_upper_bound_twice(optimizer, Float64)
end
