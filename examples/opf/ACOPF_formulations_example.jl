path = ["../examples/opf/14bus/case.raw"; "../examples/opf/14bus/case.rop"]

start_time = time();

network = process_data(path)

M = create_PowersenseData(network)

run_opf!(M, solver = Ipopt.Optimizer, obj_type = "linear", formulation = PNPAPVmodel);

@test isapprox(M.cost, 1.467e+02, rtol=1e-1)






