using ActiveSetMethods, GLPK
using JuMP

"""
The attributes to the ActiveSetMethods and the external solver can be added using the following method with JuMP
"""
#optimizer_with_attributes( ActiveSetMethods.Optimizer, "external_optimizer" => optimizer_with_attributes((GLPK.Optimizer,"msg_lev"=>GLPK.MSG_OFF)))

model = Model(ActiveSetMethods.Optimizer);
set_optimizer_attribute(model, "external_optimizer", GLPK.Optimizer)

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
JuMP.optimize!(model);

xsol = JuMP.value.(X)
ysol = JuMP.value.(Y)
status = termination_status(model)

println("Xsol = ", xsol);
println("Ysol = ", ysol);

println("Status: ", status);
