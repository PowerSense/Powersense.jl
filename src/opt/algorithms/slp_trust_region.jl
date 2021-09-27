"""
    Sequential linear programming with trust region

This is based on the paper: https://doi.org/10.1016/j.epsr.2018.09.002

# Note
- The predicted reduction ϕ_pre implemented in the original paper failed to solve many test instances. So, we have modified the predicted reduction ϕ_pre to use the directional derivative.
- We have implemented the feasibility restoration, which is not available in the original paper.
"""
mutable struct SlpTR{T,Tv,Tt} <: AbstractSlpOptimizer
    problem::Model{T,Tv,Tt} # problem data

    x::Tv # primal solution
    p::Tv # search direction
    p_slack::Dict{Int,Tv} # search direction at feasibility restoration phase
    lambda::Tv # Lagrangian dual multiplier
    mult_x_L::Tv # reduced cost for lower bound
    mult_x_U::Tv # reduced cost for upper bound

    # Evaluations at `x`
    f::T # objective function
    df::Tv # gradient
    E::Tv # constraint
    dE::Tv # Jacobian

    phi::T # merit function value
    ν::Tv # penalty parameters for the merit function

    Δ::T # current trust region size
    Δ_max::T # maximum trust region size allowed
    alpha1::T # parameter to adjust trust region size
    alpha2::T # parameter to adjust trust region size

    prim_infeas::T # primal infeasibility at `x`
    dual_infeas::T # dual (approximate?) infeasibility
    compl::T # complementary slackness

    optimizer::Union{Nothing,AbstractSubOptimizer} # Subproblem optimizer

    options::Parameters

    feasibility_restoration::Bool # indicator for feasibility restoration
    iter::Int # iteration counter
    ret::Int # solution status
    start_time::Float64 # solution start time

    function SlpTR(problem::Model{T,Tv,Tt}) where {T,Tv<:AbstractArray{T},Tt}
        slp = new{T,Tv,Tt}()
        slp.problem = problem
        slp.x = Tv(undef, problem.n)
        slp.p = zeros(problem.n)
        slp.p_slack = Dict()
        slp.lambda = zeros(problem.m)
        slp.mult_x_L = zeros(problem.n)
        slp.mult_x_U = zeros(problem.n)
        slp.df = Tv(undef, problem.n)
        slp.E = Tv(undef, problem.m)
        slp.dE = Tv(undef, length(problem.j_str))
        slp.phi = Inf

        slp.ν = Tv(undef, problem.m)
        slp.Δ = problem.parameters.tr_size
        slp.Δ_max = 2.0
        slp.alpha1 = 0.1
        slp.alpha2 = 0.25

        slp.prim_infeas = Inf
        slp.dual_infeas = Inf
        slp.compl = Inf

        slp.options = problem.parameters
        slp.optimizer = nothing

        slp.feasibility_restoration = false
        slp.iter = 1
        slp.ret = -5
        slp.start_time = 0.0
        return slp
    end
end

"""
    run!

Run the trust-region SLP algorithm
"""
function run!(slp::SlpTR)

    slp.start_time = time()

    if slp.options.OutputFlag == 1
        sparsity_val = ifelse(
            slp.problem.m > 0,
            length(slp.problem.j_str) / (slp.problem.m * slp.problem.n),
            0.0,
        )
        @printf("LP subproblem sparsity: %e\n", sparsity_val)
        add_statistics(slp.problem, "sparsity", sparsity_val)
    else
        Logging.disable_logging(Logging.Info)
    end

    # Set initial point from MOI
    @assert length(slp.x) == length(slp.problem.x)
    slp.x .= slp.problem.x
    # Adjust the initial point to satisfy the column bounds
    for i = 1:slp.problem.n
        if slp.problem.x_L[i] > -Inf
            slp.x[i] = max(slp.x[i], slp.problem.x_L[i])
        end
        if slp.problem.x_U[i] > -Inf
            slp.x[i] = min(slp.x[i], slp.problem.x_U[i])
        end
    end

    slp.iter = 1
    while true

        # evaluate function, constraints, gradient, Jacobian
        eval_functions!(slp)

        LP_time_start = time()
        # solve LP subproblem (to initialize dual multipliers)
        slp.p, slp.lambda, slp.mult_x_U, slp.mult_x_L, slp.p_slack, status =
            sub_optimize!(slp, slp.Δ)
        # @show slp.x, slp.p, slp.lambda, slp.p_slack

        add_statistics(slp.problem, "LP_time", time() - LP_time_start)

        if status ∉ [MOI.OPTIMAL, MOI.INFEASIBLE]
            @warn("Unexpected LP subproblem solution status")
            slp.ret == -3
            if norm_violations(slp, slp.x) <= slp.options.tol_infeas
                slp.ret = 6
            end
            break
        elseif status == MOI.INFEASIBLE
            if slp.feasibility_restoration == true
                @info "Failed to find a feasible direction"
                if slp.prim_infeas <= slp.options.tol_infeas
                    slp.ret = 6
                else
                    slp.ret = 2
                end
                break
            else
                slp.feasibility_restoration = true
                continue
            end
        end

        compute_nu!(slp)

        slp.prim_infeas = norm_violations(slp, Inf)
        slp.dual_infeas = KT_residuals(slp)
        slp.compl = norm_complementarity(slp)

        print_header(slp)
        print(slp)
        collect_statistics(slp)

        # Check the first-order optimality condition
        if slp.prim_infeas <= slp.options.tol_infeas &&
           slp.compl <= slp.options.tol_residual &&
           norm(slp.p, Inf) <= slp.options.tol_direction

            if slp.feasibility_restoration
                slp.feasibility_restoration = false
                slp.iter += 1
                continue
            elseif slp.dual_infeas <= slp.options.tol_residual
                slp.ret = 0
                break
            end
        end

        # Iteration counter limit
        if slp.iter >= slp.options.max_iter
            slp.ret = -1
            if slp.prim_infeas <= slp.options.tol_infeas
                slp.ret = 6
            end
            break
        end

        # step size computation
        rho = step_quality(slp)
        if slp.ret ∈ [0, 2, 6]
            break
        end
        if rho >= 0
            # update primal points
            slp.x += slp.p
        end

        slp.iter += 1
    end
    slp.problem.obj_val = slp.problem.eval_f(slp.x)
    slp.problem.status = Int(slp.ret)
    slp.problem.x .= slp.x
    slp.problem.g .= slp.E
    slp.problem.mult_g .= slp.lambda
    slp.problem.mult_x_U .= slp.mult_x_U
    slp.problem.mult_x_L .= slp.mult_x_L
    add_statistic(slp.problem, "iter", slp.iter)
end

"""
    step_quality

Evaluate the next step and adjust the trust region size
"""
function step_quality(slp::SlpTR)
    slp.phi = compute_phi(slp, slp.x, slp.p) - compute_phi(slp, slp.x)
    # @show compute_phi(slp, slp.x, slp.p), compute_phi(slp, slp.x)

    ϕ_pre = compute_derivative(slp)

    ρ = 0.0
    # @show slp.phi, ϕ_pre
    if abs(ϕ_pre) > 0.0
        ρ = slp.phi / ϕ_pre
        if ρ <= 0
            slp.Δ *= slp.alpha1
        elseif ρ <= 0.25
            slp.Δ *= slp.alpha2
        elseif ρ > 0.75
            slp.Δ = min(2 * slp.Δ, slp.Δ_max)
        end
    else
        ρ = -slp.phi
        if abs(slp.phi) < 1.e-8
            if slp.feasibility_restoration
                slp.feasibility_restoration = false
            else
                if slp.prim_infeas <= slp.options.tol_infeas
                    if slp.dual_infeas <= slp.options.tol_residual &&
                       slp.compl <= slp.options.tol_residual
                        slp.ret = 0
                    else
                        slp.ret = 6
                    end
                else
                    slp.ret = 2
                end
            end
        end
    end

    return ρ
end

"""
    compute_phi

Evaluate and return the merit function value for the current point x.
"""
compute_phi(slp::SlpTR, x::Tv) where {T,Tv<:AbstractArray{T}} =
    compute_phi(slp, x, 0.0, slp.p)

"""
    compute_phi

Evaluate and return the merit function value for the next point x + p.
"""
compute_phi(slp::SlpTR, x::Tv, p::Tv) where {T,Tv<:AbstractArray{T}} =
    compute_phi(slp, x, 1.0, p)

"""
    print_header

Print the header of iteration information.
"""
function print_header(slp::SlpTR)
    if slp.options.OutputFlag == 0
        return
    end
    if (slp.iter - 1) % 25 == 0
        @printf("  %6s", "iter")
        @printf("  %15s", "f(x_k)")
        @printf("  %15s", "Δ")
        @printf("  %14s", "|p|")
        @printf("  %14s", "inf_pr")
        @printf("  %14s", "inf_du")
        @printf("  %14s", "compl")
        @printf("  %10s", "time")
        @printf("\n")
    end
end

"""
    print

Print iteration information.
"""
function print(slp::SlpTR)
    if slp.options.OutputFlag == 0
        return
    end
    st = ifelse(slp.feasibility_restoration, "FR", "  ")
    @printf("%2s%6d", st, slp.iter)
    @printf("  %+6.8e", slp.f)
    @printf("  %+6.8e", slp.Δ)
    @printf("  %6.8e", norm(slp.p, Inf))
    @printf("  %6.8e", slp.prim_infeas)
    @printf("  %6.8e", slp.dual_infeas)
    @printf("  %6.8e", slp.compl)
    @printf("  %10.2f", time() - slp.start_time)
    @printf("\n")
end

"""
    collect_statistics

Collect iteration information.
"""
function collect_statistics(slp::SlpTR)
    if slp.options.StatisticsFlag == 0
        return
    end
    add_statistics(slp.problem, "f(x)", slp.f)
    add_statistics(slp.problem, "|p|", norm(slp.p, Inf))
    add_statistics(slp.problem, "|J|2", norm(slp.dE, 2))
    add_statistics(slp.problem, "|J|inf", norm(slp.dE, Inf))
    add_statistics(slp.problem, "inf_pr", slp.prim_infeas)
    # add_statistics(slp.problem, "inf_du", dual_infeas)
    add_statistics(slp.problem, "compl", slp.compl)
    add_statistics(slp.problem, "time_elapsed", time() - slp.start_time)
end
