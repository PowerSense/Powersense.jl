model = Model(Ipopt.Optimizer); 

path = ["../examples/opf/14bus/case.raw"; "../examples/opf/14bus/case.rop"]

start_time = time();

network = process_data(path)

M = create_PowersenseData(network, start_time)

run_OPF(M, solver = optimizer_solver, obj_type = "linear", formulation = PNPAPVmodel);

@test isapprox(M.cost, 1.467e+02, rtol=1e-1)






