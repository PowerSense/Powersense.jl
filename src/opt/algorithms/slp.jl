"""
    AbstractSlpOptimizer

Abstract type of SLP solvers
"""
abstract type AbstractSlpOptimizer <: AbstractOptimizer end

function LpData(slp::AbstractSlpOptimizer)
	A = compute_jacobian_matrix(slp)
	return QpData(
        MOI.MIN_SENSE,
        nothing,
		slp.df,
		slp.f,
		A,
		slp.E,
		slp.problem.g_L,
		slp.problem.g_U,
		slp.problem.x_L,
		slp.problem.x_U)
end

function convex_initialization!(slp::AbstractSlpOptimizer)
    if slp.options.OutputFlag == 1
        @info "Intializing with Convex Fesible Solution ..."
    end
    n = slp.problem.n
    for i = 1:n
        MOI.add_constraint(
            slp.problem.convex_model,
                    MOI.ScalarAffineFunction(
                        MOI.ScalarAffineTerm.(
                            [1.0; 1.0; -1.0],
                            [MOI.VariableIndex(i); MOI.VariableIndex(n + i); MOI.VariableIndex(2 * n + i)],
                        ),
                        0.0,
                    ),
                    MOI.EqualTo(slp.x[i]),
        )
    end
    MOI.set(slp.problem.convex_model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    
    MOI.optimize!(slp.problem.convex_model)
    for i = 1:n
        slp.x[i] = MOI.get(slp.problem.convex_model, MOI.VariablePrimal(), MOI.VariableIndex(i))
    end
    return slp.x
end

function sub_optimize!(slp::AbstractSlpOptimizer, Δ = 1000.0)
    if isnothing(slp.optimizer)
        j_row = Vector{Int}(undef, length(slp.problem.j_str))
        j_col = Vector{Int}(undef, length(slp.problem.j_str))
        for i=1:length(slp.problem.j_str)
            j_row[i] = Int(slp.problem.j_str[i][1]);
            j_col[i] = Int(slp.problem.j_str[i][2]);
        end
        slp.optimizer = QpModel(
            MOI.instantiate(slp.options.external_optimizer),
            LpData(slp),
            j_row,
            j_col,
            slp.problem.g_order
        )
        create_model!(slp.optimizer, slp.x, Δ)
    else
        slp.optimizer.data = LpData(slp)
    end
    return sub_optimize!(
        slp.optimizer,
        slp.x,
        Δ,
        slp.feasibility_restoration,
        slp.options.tol_error,
    )
end

"""
    compute_nu!

Compute the penalty parameter for the merit function. This is based on the paper: https://doi.org/10.1016/j.epsr.2018.09.002
"""
function compute_nu!(slp::AbstractSlpOptimizer)
    if slp.iter == 1
        norm_df = ifelse(slp.feasibility_restoration, 1.0, norm(slp.df))
        J = compute_jacobian_matrix(slp)
        for i = 1:slp.problem.m
            slp.ν[i] = max(1.0, norm_df / max(1.0, norm(J[i, :])))
        end
    else
        for i = 1:slp.problem.m
            slp.ν[i] = max(slp.ν[i], abs(slp.lambda[i]))
        end
    end
end

"""
    compute_phi

Evaluate and return the merit function value for a given point x + α * p.

# Arguments
- `slp`: SlpLS structure
- `x`: the current solution point
- `α`: step size taken from `x`
- `p`: direction taken from `x`
"""
function compute_phi(slp::AbstractSlpOptimizer, x::Tv, α::T, p::Tv) where {T,Tv<:AbstractArray{T}}
    ϕ = 0.0
    xp = x + α * p
    E = ifelse(α == 0.0, slp.E, slp.problem.eval_g(xp, zeros(slp.problem.m)))
    if slp.feasibility_restoration
        for i = 1:slp.problem.m
            if E[i] > slp.problem.g_U[i]
                ϕ += (E[i] - slp.problem.g_U[i])
            elseif E[i] < slp.problem.g_L[i]
                ϕ += (slp.problem.g_L[i] - E[i])
            end
        end
    else
        ϕ = slp.problem.eval_f(xp)
        for i = 1:slp.problem.m
            if E[i] > slp.problem.g_U[i]
                ϕ += slp.ν[i] * (E[i] - slp.problem.g_U[i])
            elseif E[i] < slp.problem.g_L[i]
                ϕ += slp.ν[i] * (slp.problem.g_L[i] - E[i])
            end
        end
    end
    return ϕ
end

"""
    compute_derivative

Compute the directional derivative at current solution for a given direction.
"""
function compute_derivative(slp::AbstractSlpOptimizer)
    D = 0.0
    if slp.feasibility_restoration
        for (i, v) in slp.p_slack
            D += sum(v)
        end
        for i = 1:slp.problem.m
            if slp.E[i] > slp.problem.g_U[i]
                D -= (slp.E[i] - slp.problem.g_U[i])
            elseif slp.E[i] < slp.problem.g_L[i]
                D -= (slp.problem.g_L[i] - slp.E[i])
            end
        end
    else
        D = slp.df' * slp.p
        for i = 1:slp.problem.m
            if slp.E[i] > slp.problem.g_U[i]
                D -= slp.ν[i] * (slp.E[i] - slp.problem.g_U[i])
            elseif slp.E[i] < slp.problem.g_L[i]
                D -= slp.ν[i] * (slp.problem.g_L[i] - slp.E[i])
            end
        end
    end
    return D
end

"""
    KT_residuals

Compute Kuhn-Turck residuals
"""
KT_residuals(slp::AbstractSlpOptimizer) = KT_residuals(slp.df, slp.lambda, slp.mult_x_U, slp.mult_x_L, compute_jacobian_matrix(slp))

"""
    norm_complementarity

Compute the normalized complementeraity
"""
norm_complementarity(slp::AbstractSlpOptimizer, p = Inf) = norm_complementarity(
    slp.E, slp.problem.g_L, slp.problem.g_U, 
    slp.x, slp.problem.x_L, slp.problem.x_U, 
    slp.lambda, slp.mult_x_U, slp.mult_x_L, 
    p
)

"""
    norm_violations

Compute the normalized constraint violation
"""

norm_violations(slp::AbstractSlpOptimizer, p = 1) = norm_violations(
    slp.E, slp.problem.g_L, slp.problem.g_U, 
    slp.x, slp.problem.x_L, slp.problem.x_U, 
    p
)

norm_violations(slp::AbstractSlpOptimizer, x::Tv, p = 1) where {T, Tv<:AbstractArray{T}} = norm_violations(
    slp.problem.eval_g(x, zeros(slp.problem.m)), slp.problem.g_L, slp.problem.g_U, 
    slp.x, slp.problem.x_L, slp.problem.x_U, 
    p
)

function eval_functions!(slp::AbstractSlpOptimizer)
    slp.f = slp.problem.eval_f(slp.x)
    slp.problem.eval_grad_f(slp.x, slp.df)
    slp.problem.eval_g(slp.x, slp.E)
    slp.problem.eval_jac_g(slp.x, :eval, [], [], slp.dE)
end


function evaluate_total_slack!(slp::AbstractSlpOptimizer)
    slp.slack_sum = 0.0
    for (_, v) in slp.p_slack
        slp.slack_sum += sum(v)
    end
end

compute_jacobian_matrix(slp::AbstractSlpOptimizer) = compute_jacobian_matrix(slp.problem.m, slp.problem.n, slp.problem.j_str, slp.dE)

include("slp_line_search.jl")
include("slp_trust_region.jl")
