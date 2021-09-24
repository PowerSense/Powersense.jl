using Powersense
using PowerModels, JuMP, GLPK

build_acp(data_file::String) = instantiate_model(PowerModels.parse_file(data_file), ACPPowerModel, PowerModels.build_opf)
build_acr(data_file::String) = instantiate_model(PowerModels.parse_file(data_file), ACRPowerModel, PowerModels.build_opf)
build_iv(data_file::String) = instantiate_model(PowerModels.parse_file(data_file), IVRPowerModel, PowerModels.build_opf_iv)
build_dcp(data_file::String) = instantiate_model(PowerModels.parse_file(data_file), DCPPowerModel, PowerModels.build_opf_iv)

run_opf(data_file::String) = run_opf(build_acp(data_file))
run_opf(pm::AbstractPowerModel, max_iter::Int = 100) = optimize_model!(pm, optimizer = optimizer_with_attributes(
    Powersense.Optimizer, 
    "external_optimizer" => GLPK.Optimizer,
    "max_iter" => max_iter
))

include("init_opf.jl")

#=
One can run the following:

pm = build_acp("case3.m")
init_vars(pm)
# or 
init_vars_from_ipopt(pm, build_acp("case3.m"))
run_opf(pm)

pm = build_acr("case3.m")

pm = build_iv("case3.m")

=#
