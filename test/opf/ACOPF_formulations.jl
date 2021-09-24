path = ["../examples/opf/14bus/case.raw"; "../examples/opf/14bus/case.rop"]

optimizer_solver = optimizer_with_attributes(Powersense.Optimizer,"external_optimizer"=>optimizer_with_attributes(GLPK.Optimizer,"msg_lev"=>GLPK.MSG_OFF)); 

M = create_PowersenseData(path)

run_opf!(M, solver = optimizer_solver, obj_type = "linear", formulation = PNPAPVmodel);

@test isapprox(M.cost, 1.467e+02, rtol=1e-1)






