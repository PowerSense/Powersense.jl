"""
    AbstractOptimizer

Abstract type of active set solvers
"""
abstract type AbstractOptimizer end

"""
    active_set_optimize!
    
Abstract function of running active set algorithm
"""
function active_set_optimize! end


include("algorithms/common.jl")
include("algorithms/merit.jl")
include("algorithms/subproblem.jl")


"""
    AbstractSlpOptimizer

Abstract type of SLP solvers
"""
abstract type AbstractSlpOptimizer <: AbstractOptimizer end

include("algorithms/slp_line_search.jl")
