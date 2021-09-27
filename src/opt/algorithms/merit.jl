"""
    compute_phi

Compute and return merit function value

# Arguments
- `f`: objective function value
- `mu`: penalty parameter
- `norm_E`: norm of constraint violations
"""
compute_phi(f, mu, norm_E) = f + mu * norm_E

"""
    compute_mu_merit

Compute the penalty parameter of merit function

# Arguments
- `df`: evaluation of the objective gradient
- `p`: search direction
- `rho`: parameter in (0,1)
- `norm_E`: norm of constraint violations
- `lambda`: Lagrangian multiplier
"""
function compute_mu_merit(
    df::Tv, p::Tv, rho::Float64, norm_E::T, lambda::Tv
) where {T, Tv<:AbstractArray{T}}
    mu = norm(lambda, Inf) + 1.e-4
    if norm_E > 1.e-10
        return max((df' * p) / ((1 - rho) * norm_E), mu)
    end
    return mu
end

"""
    compute_derivative

Compute and return directional derivative

# Arguments
- `df`: evaluation of the objective gradient
- `p`: search direction
- `mu`: penalty parameter
- `norm_E`: norm of constraint violations
"""
compute_derivative(df::Tv, p::Tv, mu::T, norm_E::T) where {T, Tv<:AbstractArray{T}} = df' * p - mu * norm_E
