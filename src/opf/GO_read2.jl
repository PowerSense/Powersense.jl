function process_data(path)
    if typeof(path) != String
    	(network, bus_shift)=process_raw(true, path[1]);
	network=process_rop(true, path[2], network, bus_shift);
	return network
    elseif typeof(path) == String
    	data = PowerModels.parse_file(path);
    	(network, bus_shift)=process_raw(data);
    	return network
    else
    	println("The data format is not recognized or mismatched")
    end
end

function update_format(network)
    bus_index=Dict{String,Any}();
    series = 0
    for i in keys(network["bus"]) series += 1
        bus_index=merge(bus_index, Dict(i=>series));
        network["bus"][i]=merge(network["bus"][i], Dict("adjusted_bus"=>series));
        network["bus"][i]=merge(network["bus"][i], Dict("gen"=>false));
        network["bus"][i]=merge(network["bus"][i], Dict("Pd"=>0.0));
        network["bus"][i]=merge(network["bus"][i], Dict("Qd"=>0.0));
        network["bus"][i]=merge(network["bus"][i], Dict("GL"=>0.0));
        network["bus"][i]=merge(network["bus"][i], Dict("BL"=>0.0));
        network["bus"][i]=merge(network["bus"][i], Dict("sh_max"=>0.0));
        network["bus"][i]=merge(network["bus"][i], Dict("sh_min"=>0.0));
        network["bus"][i]=merge(network["bus"][i], Dict("bcs"=>0.0));
        network["bus"][i]=merge(network["bus"][i], Dict("gfM"=>0.0));
        network["bus"][i]=merge(network["bus"][i], Dict("bfM"=>0.0));
        if network["bus"][i]["bus_type"] == 3
        	network=merge(network, Dict("slack"=>parse(Int,i)))
        end
    end
    for i in keys(network["load"])
        I_index=network["load"][i]["load_bus"];
        I=bus_index[string(I_index)];
        ST=network["load"][i]["status"];
        PL=ST*(network["load"][string(i)]["pd"])+network["bus"][string(I_index)]["Pd"];
        QL=ST*(network["load"][string(i)]["qd"])+network["bus"][string(I_index)]["Qd"];
        network["bus"][string(I_index)]=merge(network["bus"][string(I_index)], Dict("Pd"=>PL, "Qd"=>QL));
        #network["bus"][string(I)]=merge(net_bus[string(I)], Dict("Qd"=>QL));
    end
    for i in keys(network["shunt"])
        I_index=network["shunt"][i]["shunt_bus"];
        I=bus_index[string(I_index)];
        ST=network["shunt"][i]["status"];
        GL=ST*(network["shunt"][string(i)]["gs"])+network["bus"][string(I_index)]["GL"];
        BL=ST*(network["shunt"][string(i)]["bs"])+network["bus"][string(I_index)]["BL"];
        network["bus"][string(I_index)]=merge(network["bus"][string(I_index)], Dict("GL"=>GL, "BL"=>BL));
    end
    for i in keys(network["gen"])
        I_index=network["gen"][string(i)]["gen_bus"];
        I=bus_index[string(I_index)];
        network["gen"][i]=merge(network["gen"][i], Dict("gen_bus_index"=>I));
        network["gen"][i]=merge(network["gen"][i], Dict("gen_bus"=>I_index));
        network["gen"][i]=merge(network["gen"][i], Dict("gov_droop"=>0.0));
        network["gen"][i]=merge(network["gen"][i], Dict("total_apf"=>0.0));
        network["gen_cost_order"] = max(length(network["gen"][string(i)]["cost"]),network["gen_cost_order"]);
    end
    for i in keys(network["branch"])
	I_index=network["branch"][string(i)]["f_bus"];
	I=bus_index[string(I_index)];
	J_index=network["branch"][string(i)]["t_bus"];
	J=bus_index[string(J_index)];
	network["branch"][i]=merge(network["branch"][i], Dict("f_bus_index"=>I));
	network["branch"][i]=merge(network["branch"][i], Dict("t_bus_index"=>J));
	network["branch"][i]=merge(network["branch"][i], Dict("gfM"=>0.0));
	network["branch"][i]=merge(network["branch"][i], Dict("bfM"=>0.0));
    end
    return(network,bus_index)
end

function process_raw(network)
    network=merge(network, Dict("nbal"=>Dict{String,Any}(), "sshunt"=>Dict{String,Any}(), "gen_cost_order"=>0));
    network["sshunt"]=merge(network["sshunt"], Dict("bus_index"=>zeros(Int32,0)));
    
    (network, bus_shift) = update_format(network);
    
    for i=1:length(network["bus"])
        network["nbal"]=merge(network["nbal"], Dict(string(i)=>Dict()));
        network["nbal"][string(i)]=merge(network["nbal"][string(i)], Dict("gen"=>zeros(Int32,0)), Dict("Sij"=>zeros(Int32,0)), Dict("Sji"=>zeros(Int32,0)), Dict("shunt"=>zeros(Int32,0)));
    end

    for g=1:length(network["gen"]) gen=network["gen"][string(g)];
        if (gen["gen_bus_index"]*gen["gen_status"]>0) 
        push!(network["nbal"][string(gen["gen_bus_index"])]["gen"],g); 
        end
    end

    for b=1:length(network["branch"]) br=network["branch"][string(b)];
        if (br["f_bus_index"]*br["br_status"]>0) push!(network["nbal"][string(br["f_bus_index"])]["Sij"], b); end
        if (br["t_bus_index"]*br["br_status"]>0) push!(network["nbal"][string(br["t_bus_index"])]["Sji"], b); end
    end

    return(network, bus_shift)
end
