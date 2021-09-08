using JuMP
using MathOptInterface
const MOI = MathOptInterface


optimizer_solver = nothing
model = nothing
xsol = 0.0
ysol = 0.0
status = nothing

optimizer_solver = try 
			optimizer_with_attributes(ActiveSetMethods.Optimizer,"external_optimizer"=>optimizer_with_attributes(GLPK.Optimizer,"msg_lev"=>GLPK.MSG_OFF)); 
		    catch; 
		    	ActiveSetMethods.Optimizer; 
		    end

model = try 
	   Model(optimizer_solver); 
	 catch; 
	   Model(ActiveSetMethods.Optimizer); 
	 end

if typeof(optimizer_solver) == DataType 
set_optimizer_attribute(model, "external_optimizer", GLPK.Optimizer)
end

@variable(model, X);
@variable(model, Y);
@objective(model, Min, X^2 + X);
@NLconstraint(model, X^2 - X == 2);
@NLconstraint(model, X*Y == 1);
@NLconstraint(model, X*Y >= 0);
@constraint(model, X >= -2);

println("________________________________________");
print(model);
println("________________________________________");

if (typeof(optimizer_solver) == MOI.OptimizerWithAttributes) status = try 
										JuMP.optimize!(model); 
										termination_status(model) 
								         catch; 
								         	nothing 
								         end
end

if status != nothing 
xsol = JuMP.value.(X)
ysol = JuMP.value.(Y)
end






