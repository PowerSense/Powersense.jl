path = ["../examples/opf/14bus/case.raw"; "../examples/opf/14bus/case.rop"]

M = create_PowersenseData(path)

run_opf!(M, solver = Ipopt.Optimizer, obj_type = "linear", formulation = PNPAPVmodel);

@test isapprox(M.cost, 1.467e+02, rtol=1e-1)

run_opf!(M, solver = Ipopt.Optimizer, obj_type = "linear", box_constraints = true, formulation = PBRARVmodel);

@test isapprox(M.cost, 1.467e+02, rtol=1e-1)

run_opf!(M, solver = Ipopt.Optimizer, obj_type = "linear", box_constraints = false, formulation = CBRARVmodel);

@test isapprox(M.cost, 1.467e+02, rtol=1e-1)

run_opf!(M, solver = Ipopt.Optimizer, obj_type = "linear", box_constraints = true, formulation = CBRARVmodel);

@test isapprox(M.cost, 1.467e+02, rtol=1e-1)

path2 = "../examples/opf/case300.m"

M2 = create_PowersenseData(path2)

run_opf!(M2, solver = Ipopt.Optimizer, obj_type = "quadratic", formulation = PNPAPVmodel);

@test isapprox(M.cost, 7.1973e+05, rtol=1e+3)





