path = ["../examples/opf/14bus/case.raw"; "../examples/opf/14bus/case.rop"]

M = create_PowersenseData(path)

run_opf!(M, solver = Ipopt.Optimizer, obj_type = "linear", formulation = PNPAPVmodel);

@test isapprox(M.cost, 1.467e+02, rtol=1e-1)






