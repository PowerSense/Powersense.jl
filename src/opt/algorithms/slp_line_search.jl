"""
    Sequential linear programming with line search
"""
mutable struct SlpLS{T,Tv,Tt} <: AbstractSlpOptimizer
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
    directional_derivative::T # directional derivative
    alpha::T # stepsize

    prim_infeas::T # primal infeasibility at `x`
    dual_infeas::T # dual (approximate?) infeasibility
    compl::T # complementary slackness

    optimizer::Union{Nothing,AbstractSubOptimizer} # Subproblem optimizer

    options::Parameters

    feasibility_restoration::Bool # indicator for feasibility restoration
    iter::Int # iteration counter
    ret::Int # solution status
    start_time::Float64 # solution start time
    start_iter_time::Float64 # iteration start time

    function SlpLS(problem::Model{T,Tv,Tt}) where {T,Tv<:AbstractArray{T},Tt}
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

        slp.alpha = 1.0

        slp.prim_infeas = Inf
        slp.dual_infeas = Inf
        slp.compl = Inf

        slp.options = problem.parameters
        slp.optimizer = nothing

        slp.feasibility_restoration = false
        slp.iter = 1
        slp.ret = -5
        slp.start_time = 0.0
        slp.start_iter_time = 0.0

        return slp
    end
end

"""
    run!

Run the line-search SLP algorithm
"""
function run!(slp::SlpLS)

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
    is_valid_step = true
    while true

        slp.start_iter_time = time()

        # evaluate function, constraints, gradient, Jacobian
        eval_functions!(slp)
        slp.alpha = 0.0
        slp.prim_infeas = norm_violations(slp, Inf)
        slp.dual_infeas = KT_residuals(slp)
        slp.compl = norm_complementarity(slp)

        LP_time_start = time()
        # solve LP subproblem (to initialize dual multipliers)
        slp.p, slp.lambda, slp.mult_x_U, slp.mult_x_L, slp.p_slack, status =
            sub_optimize!(slp)

        add_statistics(slp.problem, "LP_time", time() - LP_time_start)

        if status ∉ [MOI.OPTIMAL, MOI.INFEASIBLE]
            @warn("Unexpected LP subproblem solution status ($status)")
            slp.ret == -3
            if slp.prim_infeas <= slp.options.tol_infeas
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
                # println("Feasibility restoration ($(status), |p| = $(norm(slp.p, Inf))) begins.")
                slp.feasibility_restoration = true
                continue
            end
        end

        compute_nu!(slp)
        slp.phi = compute_phi(slp)
        slp.directional_derivative = compute_derivative(slp)

        # step size computation
        is_valid_step = compute_alpha(slp)

        print_header(slp)
        print(slp)
        collect_statistics(slp)

        # Iteration counter limit
        if slp.iter >= slp.options.max_iter
            slp.ret = -1
            if slp.prim_infeas <= slp.options.tol_infeas
                slp.ret = 6
            end
            break
        end

        if (
            slp.prim_infeas <= slp.options.tol_infeas &&
            slp.compl <= slp.options.tol_residual
        ) || norm(slp.p, Inf) <= slp.options.tol_direction
            if slp.feasibility_restoration
                slp.feasibility_restoration = false
                slp.iter += 1
                continue
            elseif slp.dual_infeas <= slp.options.tol_residual
                slp.ret = 0
                break
            end
        end

        # Failed to find a step size
        if !is_valid_step
            @info "Failed to find a step size"
            if slp.ret == -3
                if slp.prim_infeas <= slp.options.tol_infeas
                    slp.ret = 6
                else
                    slp.ret = 2
                end
                break
            else
                slp.feasibility_restoration = true
            end

            slp.iter += 1
            continue
        end

        # update primal points
        slp.x += slp.alpha .* slp.p

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
    compute_alpha

Compute step size for line search
"""
function compute_alpha(slp::SlpLS)::Bool
    is_valid = true
    slp.alpha = 1.0
    phi_x_p = compute_phi(slp, slp.x, slp.alpha, slp.p)
    eta = slp.options.eta

    while phi_x_p > slp.phi + eta * slp.alpha * slp.directional_derivative
        # The step size can become too small.
        if slp.alpha < slp.options.min_alpha
            # @printf("Descent step cannot be computed.\n")
            # @printf("Feasibility restoration is required but not implemented yet.\n")
            if slp.feasibility_restoration
                slp.ret = -3
            end
            is_valid = false
            break
        end
        slp.alpha *= slp.options.tau
        phi_x_p = compute_phi(slp, slp.x, slp.alpha, slp.p)
    end
    # @show phi_x_p, slp.phi, slp.alpha, slp.directional_derivative, is_valid
    return is_valid
end

"""
    compute_nu!

Compute the penalty parameter for the merit function.
"""
function compute_nu!(slp::SlpLS)
    if slp.iter == 1
        for i = 1:slp.problem.m
            slp.ν[i] = abs(slp.lambda[i])
        end
    else
        for i = 1:slp.problem.m
            slp.ν[i] = max(slp.ν[i], abs(slp.lambda[i]))
        end
    end
end

"""
    compute_phi

Evaluate and return the merit function value at the current point x
"""
compute_phi(slp::SlpLS) = compute_phi(slp, slp.x, 0.0, slp.p)

"""
    print_header

Print the header of iteration information.
"""
function print_header(slp::SlpLS)
    if slp.options.OutputFlag == 0
        return
    end
    if (slp.iter - 1) % 25 == 0
        @printf("  %6s", "iter")
        @printf("  %15s", "f(x_k)")
        # @printf("  %15s", "ϕ(x_k)")
        @printf("  %15s", "D(ϕ,p)")
        # @printf("  %15s", "∇f^Tp")
        @printf("  %14s", "α")
        @printf("  %14s", "|p|")
        @printf("  %14s", "α|p|")
        # @printf("  %14s", "|∇f|")
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
function print(slp::SlpLS)
    if slp.options.OutputFlag == 0
        return
    end
    st = ifelse(slp.feasibility_restoration, "FR", "  ")
    @printf("%2s%6d", st, slp.iter)
    @printf("  %+6.8e", slp.f)
    # @printf("  %+6.8e", slp.phi)
    @printf("  %+.8e", slp.directional_derivative)
    # @printf("  %+.8e", slp.df' * slp.p)
    @printf("  %6.8e", slp.alpha)
    @printf("  %6.8e", norm(slp.p, Inf))
    @printf("  %6.8e", slp.alpha * norm(slp.p, Inf))
    # @printf("  %.8e", norm(slp.df))
    @printf("  %6.8e", slp.prim_infeas)
    @printf("  %.8e", slp.dual_infeas)
    @printf("  %6.8e", slp.compl)
    @printf("  %10.2f", time() - slp.start_time)
    @printf("\n")
end

"""
    collect_statistics

Collect iteration information.
"""
function collect_statistics(slp::SlpLS)
    if slp.options.StatisticsFlag == 0
        return
    end
    mode = (slp.feasibility_restoration) ? "FR" : "OPT"
    add_statistics(slp.problem, "f(x)", slp.f)
    add_statistics(slp.problem, "ϕ(x_k))", slp.phi)
    add_statistics(slp.problem, "D(ϕ,p)", slp.directional_derivative)
    add_statistics(slp.problem, "|p|", norm(slp.p, Inf))
    add_statistics(slp.problem, "|J|2", norm(slp.dE, 2))
    add_statistics(slp.problem, "|J|inf", norm(slp.dE, Inf))
    add_statistics(slp.problem, "inf_pr", slp.prim_infeas)
    add_statistics(slp.problem, "mode", mode)
    # add_statistics(slp.problem, "inf_du", dual_infeas)
    add_statistics(slp.problem, "compl", slp.compl)
    add_statistics(slp.problem, "alpha", slp.alpha)
    add_statistics(slp.problem, "iter_time", time() - slp.start_iter_time)
    add_statistics(slp.problem, "time_elapsed", time() - slp.start_time)
end