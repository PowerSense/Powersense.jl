struct QpData{T, Tv<:AbstractArray{T}, Tm<:AbstractMatrix{T}}
	sense::MOI.OptimizationSense
	#Q::Union{Nothing,Tm}
	Q::Union{Nothing,Tm}
	c::Tv
	c0::T # objective functiton constant term
	A::Tm
	b::Tv
	c_lb::Tv
	c_ub::Tv
	v_lb::Tv
	v_ub::Tv
	j_row::Array{Int64,1}
    	j_col::Array{Int64,1}
    	adj::Array{Int64,1}
	constr_v_ub
	constr_v_lb
	constr
end

"""
	sub_optimize!

Solve subproblem

# Arguments
- `model`: MOI abstract optimizer
- `qp`: QP problem data
- `mu`: penalty parameter
- `x_k`: trust region center
- `Δ`: trust region size
"""

function sub_optimize!(
	model::MOI.AbstractOptimizer,
	qp::QpData{T,Tv,Tm},
	mu::T,
	x_k::Tv,
	Δ::T,
	tol_error::T,
	condition_flag::Int
) where {T, Tv, Tm}
	
	qp.c[abs.(qp.c) .<= tol_error] .= zero(eltype(qp.c))
	
	n = length(qp.c);
	m = length(qp.c_lb);
	
	@assert n > 0
	@assert m >= 0
	@assert length(qp.c) == n
	@assert length(qp.c_lb) == m
	@assert length(qp.c_ub) == m
	@assert length(qp.v_lb) == n
	@assert length(qp.v_ub) == n
	@assert length(x_k) == n
	
	MOI.modify(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), MOI.ScalarConstantChange(qp.c0))
	for i in 1:n
		MOI.modify(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), MOI.ScalarCoefficientChange(MOI.VariableIndex(i), qp.c[i]))
	end


	# Slacks v and u are added only for constrained problems.
	if m > 0
		for i in 1:m
			MOI.modify(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), MOI.ScalarCoefficientChange(MOI.VariableIndex(i+n), mu))
			MOI.modify(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), MOI.ScalarCoefficientChange(MOI.VariableIndex(i+n+m), mu))
		end
	end

	# set constant term to the objective function
	for i = 1:n
		ub = min(Δ, qp.v_ub[i] - x_k[i])
		lb = max(-Δ, qp.v_lb[i] - x_k[i])
		ub = (abs(ub) <= tol_error) ? 0.0 : ub;
		lb = (abs(lb) <= tol_error) ? 0.0 : lb;
		MOI.set(model, MOI.ConstraintSet(), qp.constr_v_ub[i], MOI.LessThan(ub))
		MOI.set(model, MOI.ConstraintSet(), qp.constr_v_lb[i], MOI.GreaterThan(lb))
	end

	#condition1 = cond(Array(qp.A),1)
	#condition = cond(Array(qp.A))
	#conditionInf = cond(Array(qp.A),Inf)
	#@show condition
	for i=1:length(qp.j_row)
		coeff = (abs(qp.A[qp.j_row[i],qp.j_col[i]]) <= tol_error) ? 0.0 : qp.A[qp.j_row[i],qp.j_col[i]];
		MOI.modify(model, qp.constr[qp.j_row[i]], MOI.ScalarCoefficientChange(MOI.VariableIndex(qp.j_col[i]), coeff))	
	end
	
	for ind in 1:length(qp.adj)
		i_row = qp.adj[ind]
		for i in 1:length(qp.A[i_row,:].nzind)
			i_col = qp.A[i_row,:].nzind[i];
			value = qp.A[i_row,:].nzval[i];
			coeff = (abs(value) <= tol_error) ? 0.0 : value;
			MOI.modify(model, qp.constr[m+ind], MOI.ScalarCoefficientChange(MOI.VariableIndex(i_col), coeff));
		end		
	end
	
	for i=1:m
		c_ub = qp.c_ub[i]-qp.b[i];
		c_lb = qp.c_lb[i]-qp.b[i];
		
		c_ub = (abs(c_ub) <= tol_error) ? 0.0 : c_ub;
		c_lb = (abs(c_lb) <= tol_error) ? 0.0 : c_lb;
		
		
		if qp.c_lb[i] == qp.c_ub[i] #This means the constraint is equality
			MOI.set(model, MOI.ConstraintSet(), qp.constr[i], MOI.EqualTo(c_lb))
		elseif qp.c_lb[i] != -Inf && qp.c_ub[i] != Inf && qp.c_lb[i] < qp.c_ub[i]
			MOI.set(model, MOI.ConstraintSet(), qp.constr[i], MOI.GreaterThan(c_lb))
		elseif qp.c_lb[i] != -Inf
			MOI.set(model, MOI.ConstraintSet(), qp.constr[i], MOI.GreaterThan(c_lb))
		elseif qp.c_ub[i] != Inf
			MOI.set(model, MOI.ConstraintSet(), qp.constr[i], MOI.LessThan(c_ub))
		end
	end
	
	for i in 1:length(qp.adj) 
		c_ub = qp.c_ub[qp.adj[i]]-qp.b[qp.adj[i]];	
		c_ub = (abs(c_ub) <= tol_error) ? 0.0 : c_ub;
		MOI.set(model, MOI.ConstraintSet(), qp.constr[i+m], MOI.LessThan(c_ub))
	end

	MOI.optimize!(model)
	#println("Condition number: ", MOI.get(model, Gurobi.ModelAttribute("KappaExact")))
	TerminationStatus = MOI.get(model, MOI.TerminationStatus())
	PrimalStatus = MOI.get(model, MOI.PrimalStatus())
	ResultCount = MOI.get(model, MOI.ResultCount())

	# TODO: These can be part of data.
	Xsol = Tv(undef, n)
	if m > 0
		Usol = Tv(undef, m);
		Vsol = Tv(undef, m);
	end
	lambda = Tv(undef, m)
	mult_x_U = Tv(undef, n)
	mult_x_L = Tv(undef, n)
	infeasibility = 0.0
   	
	if PrimalStatus == MOI.FEASIBLE_POINT
		for i=1:n
			Xsol[i] = MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(i));
		end
		if m > 0
			for i=1:m
				Usol[i] = MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(i+n));
				Vsol[i] = MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(i+n+m));
			end
		end

		# extract the multipliers to constraints
		lambda_indv = MOI.get(model, MOI.ConstraintDual(1), qp.constr)
		lambda .= lambda_indv[1:m];
		for ind in 1:length(qp.adj)
			i_row = qp.adj[ind]
			lambda[i_row] += lambda_indv[m+ind];
		end

		# extract the multipliers to column bounds
		mult_x_U .= MOI.get(model, MOI.ConstraintDual(1), qp.constr_v_ub)
		mult_x_L .= MOI.get(model, MOI.ConstraintDual(1), qp.constr_v_lb)
		# careful because of the trust region
		for j=1:n
			if Xsol[j] < qp.v_ub[j] - x_k[j]
				mult_x_U[j] = 0.0
			end
			if Xsol[j] > qp.v_lb[j] - x_k[j]
				mult_x_L[j] = 0.0
			end
		end
	elseif TerminationStatus == MOI.DUAL_INFEASIBLE
		@error "Trust region must be employed."
	else
		@error "Unexpected status: $(TerminationStatus)"
	end
	
	simplex = try 
			MOI.get(model, MOI.SimplexIterations()); 
		    catch; 
		    	0.0; 
		    end
		    
	barrier = try 
			MOI.get(model, MOI.BarrierIterations()); 
		    catch; 
		    	0.0; 
		    end
	
	condition = 0.0;
	
	if condition_flag > 0
	condition = try 
			MOI.get(model, Gurobi.ModelAttribute("KappaExact"));
		    catch; 
		    	0.0; 
		    end
	end
	return Xsol, lambda, mult_x_U, mult_x_L, infeasibility, PrimalStatus, simplex, barrier, condition
end


function create_model!(
	model::MOI.AbstractOptimizer,
	qp::QpData{T,Tv,Tm},
	mu::T,
	x_k::Tv,
	Δ::T,
	tol_error::T
) where {T, Tv,Tm}

	# empty optimizer just in case
	MOI.empty!(model)
	
	n = length(qp.c);
	m = length(qp.c_lb);

	@assert n > 0
	@assert m >= 0
	@assert length(qp.c) == n
	@assert length(qp.c_lb) == m
	@assert length(qp.c_ub) == m
	@assert length(qp.v_lb) == n
	@assert length(qp.v_ub) == n
	@assert length(x_k) == n
	
	# variables
	x = MOI.add_variables(model, n)
	
	# objective function
	obj_terms = Array{MOI.ScalarAffineTerm{T},1}();
	for i in 1:n
		push!(obj_terms, MOI.ScalarAffineTerm{T}(qp.c[i], MOI.VariableIndex(i)));
	end

	# Slacks v and u are added only for constrained problems.
	if m > 0
		u = MOI.add_variables(model, m)
		v = MOI.add_variables(model, m)
		append!(obj_terms, MOI.ScalarAffineTerm.(mu, u));
		append!(obj_terms, MOI.ScalarAffineTerm.(mu, v));
	end

	# set constant term to the objective function
	MOI.set(model,
		MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
		MOI.ScalarAffineFunction(obj_terms, qp.c0))
	MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)


	for i = 1:n
		ub = min(Δ, qp.v_ub[i] - x_k[i])
		lb = max(-Δ, qp.v_lb[i] - x_k[i])
		ub = (abs(ub) <= tol_error) ? 0.0 : ub;
		lb = (abs(lb) <= tol_error) ? 0.0 : lb;
		push!(qp.constr_v_ub, MOI.add_constraint(model, MOI.SingleVariable(x[i]), MOI.LessThan(ub)))
		push!(qp.constr_v_lb, MOI.add_constraint(model, MOI.SingleVariable(x[i]), MOI.GreaterThan(lb)))
	end
	
	for i=1:m
		terms = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0;-1.0], [u[i];v[i]]), 0.0);
		c_ub = qp.c_ub[i]-qp.b[i];
		c_lb = qp.c_lb[i]-qp.b[i];
		
		c_ub = (abs(c_ub) <= tol_error) ? 0.0 : c_ub;
		c_lb = (abs(c_lb) <= tol_error) ? 0.0 : c_lb;
		
		if qp.c_lb[i] == qp.c_ub[i] #This means the constraint is equality
			push!(qp.constr, MOI.add_constraint(model, terms, MOI.EqualTo(c_lb)))

		elseif qp.c_lb[i] != -Inf && qp.c_ub[i] != Inf && qp.c_lb[i] < qp.c_ub[i]
			push!(qp.constr, MOI.add_constraint(model, terms, MOI.GreaterThan(c_lb)))
			push!(qp.adj,i);
		elseif qp.c_lb[i] != -Inf
			push!(qp.constr, MOI.add_constraint(model, terms, MOI.GreaterThan(c_lb)))
		elseif qp.c_ub[i] != Inf
			push!(qp.constr, MOI.add_constraint(model, terms, MOI.LessThan(c_ub)))
		end
		MOI.add_constraint(model, MOI.SingleVariable(u[i]), MOI.GreaterThan(0.0))
		MOI.add_constraint(model, MOI.SingleVariable(v[i]), MOI.GreaterThan(0.0))
	end
	
	for i in qp.adj
		terms = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0;-1.0], [u[i];v[i]]), 0.0);
		c_ub = qp.c_ub[i]-qp.b[i];	
		c_ub = (abs(c_ub) <= tol_error) ? 0.0 : c_ub;
		push!(qp.constr, MOI.add_constraint(model, terms, MOI.LessThan(c_ub)))
	end
   return qp
end

"""
	get_moi_constraint_row_terms

Get the array of MOI constraint terms from row `i` of matrix `A`
"""
function get_moi_constraint_row_terms(A::Tm, i::Int) where {T, Tm<:AbstractSparseMatrix{T,Int}}
	Ai = A[i,:]
	terms = Array{MOI.ScalarAffineTerm{T},1}()
	for (ind, val) in zip(Ai.nzind, Ai.nzval)
		push!(terms, MOI.ScalarAffineTerm{T}(val, MOI.VariableIndex(ind)))
	end
	return terms
end

function get_moi_constraint_row_terms(A::Tm, i::Int) where {T, Tm<:DenseArray{T,2}}
	terms = Array{MOI.ScalarAffineTerm{T},1}()
	for j in 1:size(A,2)
		if !isapprox(A[i,j], 0.0)
			push!(terms, MOI.ScalarAffineTerm{T}(A[i,j], MOI.VariableIndex(j)))
		end
	end
	return terms
end
