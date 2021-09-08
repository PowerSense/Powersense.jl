using Printf

mutable struct SlpLS{T,Tv,Tt} <: AbstractSlpOptimizer
    problem::Model{T,Tv,Tt}

    x::Tv # primal solution
    p::Tv
    lambda::Tv
    mult_x_L::Tv
    mult_x_U::Tv
    
    f::T
    df::Tv
    E::Tv
    dE::Tv
    # hLag::Tv
    phi::T
    directional_derivative::T
    
    j_row::Array{Int64,1}
    j_col::Array{Int64,1}

    norm_E::T # norm of constraint violations
    mu_merit::T
    mu_lp::T
    alpha::T

    optimizer::MOI.AbstractOptimizer

    options::Parameters
    statistics::Dict{String,Any}

    iter::Int
    ret::Int

    function SlpLS(problem::Model{T,Tv,Tt}) where {T, Tv<:AbstractArray{T}, Tt}
        slp = new{T,Tv,Tt}()
        slp.problem = problem
        slp.x = Tv(undef, problem.n)
        slp.p = zeros(problem.n)
        slp.lambda = zeros(problem.m)
        slp.mult_x_L = zeros(problem.n)
        slp.mult_x_U = zeros(problem.n)
        slp.f = 0.0
        slp.df = Tv(undef, problem.n)
        slp.E = Tv(undef, problem.m)
        slp.dE = Tv(undef, length(problem.j_str))
        # slp.hLag = Tv(undef, length(problem.h_str))
        slp.phi = Inf
        
        slp.j_row = zeros(Int,length(slp.problem.j_str));
        slp.j_col = zeros(Int,length(slp.problem.j_str));
        for i=1:length(slp.problem.j_str)
             slp.j_row[i] = Int(slp.problem.j_str[i][1]);
             slp.j_col[i] = Int(slp.problem.j_str[i][2]);
    	end

        slp.norm_E = 0.0
        slp.mu_merit = problem.parameters.mu_merit
        slp.mu_lp = problem.parameters.mu_lp
        slp.alpha = 1.0

        slp.options = problem.parameters
        slp.optimizer = MOI.instantiate(slp.options.external_optimizer)

        slp.iter = 1
        slp.ret = -5
        return slp
    end
end

function active_set_optimize!(slp::SlpLS)
    Δ = slp.options.tr_size

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
    if slp.options.StatisticsFlag != 0
    	slp.problem.statistics["f(x)"] = Array{Float64,1}()
    	slp.problem.statistics["ϕ(x_k))"] = Array{Float64,1}()
    	slp.problem.statistics["D(ϕ,p)"] = Array{Float64,1}()
    	slp.problem.statistics["|p|"] = Array{Float64,1}()
    	slp.problem.statistics["|J|2"] = Array{Float64,1}() 
    	slp.problem.statistics["|J|inf"] = Array{Float64,1}() 
    	slp.problem.statistics["iter_time"] = Array{Float64,1}() 
    	slp.problem.statistics["alpha"] = Array{Float64,1}()
    	slp.problem.statistics["μ_merit"] = Array{Float64,1}()
    	slp.problem.statistics["μ_lp"] = Array{Float64,1}()
    	slp.problem.statistics["m(p)"] = Array{Float64,1}()
    	slp.problem.statistics["inf_pr"] = Array{Float64,1}()
    	slp.problem.statistics["inf_du"] = Array{Float64,1}()
    	slp.problem.statistics["compl"] = Array{Float64,1}()
    	slp.problem.statistics["Sparsity"] = Array{Float64,1}()
    	slp.problem.statistics["LP_time"] = Array{Float64,1}()
    	slp.problem.statistics["iter_time"] = Array{Float64,1}()    
    	slp.problem.statistics["LP_simplex_iter"] = Array{Float64,1}()
    	slp.problem.statistics["LP_barrier_iter"] = Array{Float64,1}()
    	slp.problem.statistics["condition"] = Array{Float64,1}()	
    end

    eval_functions!(slp)
    qp = create_model!(slp, Δ)
    
    while true
	iter_time_start = time();
        if (slp.iter - 1) % 100 == 0 && slp.options.OutputFlag != 0
            @printf("%6s", "iter")
            @printf("  %15s", "f(x_k)")
            @printf("  %15s", "ϕ(x_k)")
            @printf("  %15s", "D(ϕ,p)")
            @printf("  %15s", "∇f^Tp")
            @printf("  %14s", "|p|")
            # @printf("  %14s", "|∇f|")
            @printf("  %14s", "μ_merit")
            @printf("  %14s", "μ_lp")
            @printf("  %14s", "m(p)")
            @printf("  %14s", "inf_pr")
            @printf("  %14s", "inf_du")
            @printf("  %14s", "compl")
            @printf("  %14s", "Sparsity")
            @printf("\n")
        end
        eval_functions!(slp)
        if slp.iter != 1
        	eval_functions!(slp)
        end
        slp.norm_E = norm_violations(slp)
        LP_time_start = time()
        # solve LP subproblem (to initialize dual multipliers)
        slp.p, lambda, mult_x_U, mult_x_L, infeasibility, status, LP_simplex_iter, LP_barrier_iter, condition = sub_optimize!(slp, Δ, qp)
        # @show slp.lambda
        if slp.options.StatisticsFlag != 0
	    	push!(slp.problem.statistics["LP_time"],time()-LP_time_start)
    	end
    	
        # update multipliers
        slp.lambda .= lambda
        slp.mult_x_U .= mult_x_U
        slp.mult_x_L .= mult_x_L

        compute_mu_merit!(slp)
        slp.phi = compute_phi(slp)
        slp.directional_derivative = compute_derivative(slp)
        prim_infeas = norm_violations(slp, Inf)
        
        dual_infeas = KT_residuals(slp)
        compl = norm_complementarity(slp)
        sparsity_val = slp.problem.m > 0 ? length(slp.problem.j_str) / (slp.problem.m * slp.problem.n) : 0.0
        
        if slp.options.OutputFlag != 0
		@printf("%6d", slp.iter)
		@printf("  %+.8e", slp.f)
		@printf("  %+.8e", slp.phi)
		@printf("  %+.8e", slp.directional_derivative)
		@printf("  %+.8e", slp.df' * slp.p)
		@printf("  %.8e", norm(slp.p))
		# @printf("  %.8e", norm(slp.df))
		@printf("  %.8e", slp.mu_merit)
		@printf("  %.8e", slp.mu_lp)
		@printf("  %.8e", infeasibility)
		@printf("  %.8e", prim_infeas)
		@printf("  %.8e", dual_infeas)
		@printf("  %.8e", compl)
		@printf("  %.8e", sparsity_val)
		@printf("\n")
        end
        
        if slp.options.StatisticsFlag != 0
	    	push!(slp.problem.statistics["f(x)"],slp.f)
	    	push!(slp.problem.statistics["ϕ(x_k))"],slp.phi)
	    	push!(slp.problem.statistics["D(ϕ,p)"],slp.directional_derivative)
	    	push!(slp.problem.statistics["|p|"],norm(slp.p))
	    	push!(slp.problem.statistics["|J|2"],norm(slp.dE, 2))
	    	push!(slp.problem.statistics["|J|inf"],norm(slp.dE, Inf))
	    	push!(slp.problem.statistics["μ_merit"],slp.mu_merit)
	    	push!(slp.problem.statistics["μ_lp"],slp.mu_lp)
	    	push!(slp.problem.statistics["m(p)"],infeasibility)
	    	push!(slp.problem.statistics["inf_pr"],prim_infeas)
	    	push!(slp.problem.statistics["inf_du"],dual_infeas)
	    	push!(slp.problem.statistics["compl"],compl)
	    	push!(slp.problem.statistics["Sparsity"],sparsity_val)
	    	push!(slp.problem.statistics["LP_simplex_iter"],LP_simplex_iter)
	    	push!(slp.problem.statistics["LP_barrier_iter"],LP_barrier_iter)
	    	push!(slp.problem.statistics["condition"],condition)
    	end

        # LP subproblem is infeasible or Numerical issues
        if status != MOI.FEASIBLE_POINT
            slp.ret = 2
            break
        end
        
        # If the LP subproblem is infeasible, increase mu_lp and resolve.
        if infeasibility > 1.e-10 && slp.mu_lp < slp.options.max_mu
            slp.mu_lp = max(slp.options.max_mu, slp.mu_lp*10)
            slp.iter += 1
            continue
        end

        # Check the first-order optimality condition
        if prim_infeas <= slp.options.tol_infeas &&
            dual_infeas <= slp.options.tol_residual &&
            compl <= slp.options.tol_residual
            slp.ret = 0
            break
        end
        
        # Iteration counter limit
        if slp.iter >= slp.options.max_iter
            slp.ret = -1
            if norm_violations(slp, slp.x) <= slp.options.tol_infeas
                slp.ret = 6
            end
            break
        end
        

        # step size computation
        is_valid_step = compute_alpha(slp)
        if slp.ret == -3
            @warn "Failed to find a step size"
            break
        end
        if !is_valid_step && slp.iter < slp.options.max_iter
            slp.iter += 1
            if slp.options.StatisticsFlag != 0
	    	push!(slp.problem.statistics["alpha"],0.0)
    	    end
            continue
        end
        if slp.options.StatisticsFlag != 0
	    	push!(slp.problem.statistics["alpha"],slp.alpha)
    	end
        #@show slp.alpha

        # update primal points
        slp.x += slp.alpha .* slp.p
        
        if slp.options.StatisticsFlag != 0
	    push!(slp.problem.statistics["iter_time"],time()-iter_time_start);
    	end
        
        slp.iter += 1
    end
    if slp.options.StatisticsFlag != 0
	slp.problem.statistics["iter"] = slp.iter
    end
    #@show slp.lambda
    slp.problem.obj_val = slp.problem.eval_f(slp.x)
    slp.problem.status = Int(slp.ret)
    slp.problem.x .= slp.x
    slp.problem.g .= slp.E
    slp.problem.mult_g .= slp.lambda
    slp.problem.mult_x_U .= slp.mult_x_U
    slp.problem.mult_x_L .= slp.mult_x_L
end

function LpData(slp::SlpLS, qp=nothing)
	if qp == nothing
		constr_v_ub = MOI.ConstraintIndex[];
		constr_v_lb = MOI.ConstraintIndex[];
		constr = MOI.ConstraintIndex[];
		adj = zeros(Int,0);
	else
		constr_v_ub = qp.constr_v_ub;
		constr_v_lb = qp.constr_v_lb;
		constr = qp.constr
		adj = qp.adj
	end
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
		slp.problem.x_U,
		slp.j_row,
		slp.j_col,
		adj,
		constr_v_ub,
		constr_v_lb,
		constr)
end


sub_optimize!(slp::SlpLS, Δ, qp) = sub_optimize!(
	slp.optimizer,
	LpData(slp, qp),
	slp.mu_lp,
	slp.x,
	Δ,
	slp.options.tol_error,
	slp.options.condition_flag
)

create_model!(slp::SlpLS, Δ) = create_model!(
	slp.optimizer,
	LpData(slp),
	slp.mu_lp,
	slp.x,
	Δ,
	slp.options.tol_error
)

"""
    KT_residuals

Compute Kuhn-Turck residuals
"""
KT_residuals(slp::SlpLS) = KT_residuals(slp.df, slp.lambda, slp.mult_x_U, slp.mult_x_L, compute_jacobian_matrix(slp))

"""
    norm_complementarity

Compute the normalized complementeraity
"""
norm_complementarity(slp::SlpLS, p = Inf) = norm_complementarity(
    slp.E, slp.problem.g_L, slp.problem.g_U, 
    slp.x, slp.problem.x_L, slp.problem.x_U, 
    slp.lambda, slp.mult_x_U, slp.mult_x_L, 
    p
)

"""
    norm_violations

Compute the normalized constraint violation
"""

norm_violations(slp::SlpLS, p = 1) = norm_violations(
    slp.E, slp.problem.g_L, slp.problem.g_U, 
    slp.x, slp.problem.x_L, slp.problem.x_U, 
    p
)

norm_violations(slp::SlpLS{T,Tv,Tt}, x::Tv, p = 1) where {T, Tv<:AbstractArray{T}, Tt} = norm_violations(
    slp.problem.eval_g(x, zeros(slp.problem.m)), slp.problem.g_L, slp.problem.g_U, 
    slp.x, slp.problem.x_L, slp.problem.x_U, 
    p
)

function eval_functions!(slp::SlpLS)
    slp.f = slp.problem.eval_f(slp.x)
    slp.problem.eval_grad_f(slp.x, slp.df)
    slp.problem.eval_g(slp.x, slp.E)
    slp.problem.eval_jac_g(slp.x, :eval, [], [], slp.dE)
    # obj_factor = 1.0
    # slp.problem.eval_h(slp.x, :eval, [], [], obj_factor, slp.lambda, slp.hLag)
    # @show slp.f, slp.df, slp.E, slp.dE
end

#compute_jacobian_matrix(slp::SlpLS) = compute_jacobian_matrix(slp.problem.m, slp.problem.n, slp.problem.j_str, slp.dE)
compute_jacobian_matrix(slp::SlpLS) = compute_jacobian_matrix(slp.problem.m, slp.problem.n, slp.j_row, slp.j_col, slp.dE)

compute_mu_merit(slp::SlpLS) = compute_mu_merit(slp.df, slp.p, slp.options.rho, slp.norm_E, slp.lambda)
function compute_mu_merit!(slp::SlpLS)
    slp.mu_merit = max(slp.mu_merit, compute_mu_merit(slp))
end

"""
    compute_alpha

Compute step size for line search
"""
function compute_alpha(slp::SlpLS)::Bool
    is_valid = true

    slp.alpha = 1.0
    phi_x_p = compute_phi(slp, slp.x .+ slp.alpha * slp.p)
    eta = slp.options.eta

    while phi_x_p > slp.phi + eta * slp.alpha * slp.directional_derivative
        # The step size can become too small.
        if slp.alpha < slp.options.min_alpha
            # This can happen when mu is not sufficiently large.
            # A descent step is valid only for a sufficiently large mu.
            # slp.mu_merit = slp.mu_merit * 1.5
            # is_valid = false
            
            if slp.mu_merit < slp.options.max_mu
                # Increase mu, if mu is not large enough.
                slp.mu_merit = min(slp.options.max_mu, slp.mu_merit * 10)
                @printf("* step size too small: increase mu to %e\n", slp.mu_merit)
                is_valid = false
            # elseif eta > 1.e-6
            #     eta *= 0.5
            #     @printf("* step size too small: decrease eta to %e\n", eta)
            #     continue
            else
                slp.ret = -3
            end
            
            break
        end
        slp.alpha *= slp.options.tau
        phi_x_p = compute_phi(slp, slp.x + slp.alpha * slp.p)
        # @show phi_x_p, slp.phi, slp.alpha, slp.directional_derivative, slp.phi + slp.options.eta * slp.alpha * slp.directional_derivative
    end
    return is_valid
end

# merit function
compute_phi(slp::SlpLS) = compute_phi(slp.f, slp.mu_merit, slp.norm_E)
compute_phi(slp::SlpLS{T,Tv,Tt}, x::Tv) where {T, Tv<:AbstractArray{T}, Tt} = compute_phi(
    slp.problem.eval_f(x), slp.mu_merit, norm_violations(slp, x))

# directional derivative
compute_derivative(slp::SlpLS) = compute_derivative(slp.df, slp.p, slp.mu_merit, slp.norm_E)
