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

function convex_optimize!(slp::AbstractSlpOptimizer)
    if slp.options.OutputFlag == 1
        @info "Intializing with Convex Fesible Solution ..."
    end
    model = MOI.instantiate(slp.options.external_optimizer)
    # @show slp.x
    MOI.copy_to(model, slp.problem.convex_model)
    n = slp.problem.n
    @show n
    m = slp.problem.m
    @show m
    u = MOI.add_variables(model, n)
    v = MOI.add_variables(model, n)
    for i = 1:n
        MOI.modify(
            model,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ScalarCoefficientChange(MOI.VariableIndex(n + i), 1.0),
        )
        MOI.add_constraint(model, MOI.SingleVariable(MOI.VariableIndex(n + i)), MOI.GreaterThan(0.0))
        MOI.modify(
            model,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ScalarCoefficientChange(MOI.VariableIndex(2 * n + i), 1.0),
        )
        MOI.add_constraint(model, MOI.SingleVariable(MOI.VariableIndex(2 * n + i)), MOI.GreaterThan(0.0))
        MOI.add_constraint(
                    model,
                    MOI.ScalarAffineFunction(
                        MOI.ScalarAffineTerm.(
                            [-1.0; 1.0; -1.0],
                            [MOI.VariableIndex(i); MOI.VariableIndex(n + i); MOI.VariableIndex(2 * n + i)],
                        ),
                        slp.x[i],
                    ),
                    MOI.EqualTo(0.0),
        )
    end
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    
    MOI.optimize!(model)
    for i = 1:n
        slp.x[i] = MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(i))
    end
    # @show slp.x
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
            j_col
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
        p_slack = slp.p_slack
        ϕ = slp.prim_infeas
        for (i, v) in p_slack
            ϕ += α * sum(v)
        end
        for i = 1:slp.problem.m
            viol = maximum([0.0, slp.E[i] - slp.problem.g_U[i], slp.problem.g_L[i] - slp.E[i]])
            lhs = E[i] - viol
            if slp.problem.g_L[i] > -Inf && slp.problem.g_U[i] < Inf
                lhs += α * (p_slack[i][1] - p_slack[i][2])
            elseif slp.problem.g_L[i] > -Inf
                lhs += α * p_slack[i][1]
            elseif slp.problem.g_U[i] < Inf
                lhs -= α * p_slack[i][1]
            end
            ϕ +=
                slp.ν[i] *
                maximum([0.0, lhs - slp.problem.g_U[i], slp.problem.g_L[i] - lhs])
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
            viol = maximum([0.0, slp.E[i] - slp.problem.g_U[i], slp.problem.g_L[i] - slp.E[i]])
            lhs = slp.E[i] - viol
            D -=
                slp.ν[i] *
                maximum([0.0, lhs - slp.problem.g_U[i], slp.problem.g_L[i] - lhs])
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

compute_jacobian_matrix(slp::AbstractSlpOptimizer) = compute_jacobian_matrix(slp.problem.m, slp.problem.n, slp.problem.j_str, slp.dE)

include("slp_line_search.jl")
include("slp_trust_region.jl")
