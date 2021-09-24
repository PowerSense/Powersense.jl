include("GO_read.jl");
include("GO_read2.jl");

mutable struct PowersenseData
    nbus::Int  # Number of buses
    ngen::Int  # Number of generators
    nbr::Int  # Number of branches
    nbss::Int  # Number of controllable shunt
    Vmax::Vector{Float64}
    Vmin::Vector{Float64}
    Pmax::Vector{Float64}
    Pmin::Vector{Float64}
    Qmax::Vector{Float64}
    Qmin::Vector{Float64}
    bmax::Vector{Float64}
    bmin::Vector{Float64}
    Imax::Vector{Float64}

    Pg::Vector{Float64}
    Qg::Vector{Float64}
    Bss::Vector{Float64}
    V::Vector{Float64}
    θ::Vector{Float64}
    Gf::Vector{Float64}
    Bf::Vector{Float64}
    Gt::Vector{Float64}
    Bt::Vector{Float64}
    gf::Vector{Float64}
    bf::Vector{Float64}
    gt::Vector{Float64}
    bt::Vector{Float64}
    

    vr::Vector{Float64}
    vi::Vector{Float64}
    br::Vector{Tuple{Int64,Int64}}
    node::Dict{String,Any}
    Pd::Vector{Float64}
    Qd::Vector{Float64}
    Gl::Vector{Float64}
    Bl::Vector{Float64}

    P::Vector{Float64}
    Q::Vector{Float64}
    Iij::Vector{Float64}
    Iji::Vector{Float64}
    Pij::Vector{Float64}
    Qij::Vector{Float64}
    lps::Vector{Int}
    pgh::Vector{Vector{Float64}}
    cgh::Vector{Vector{Float64}}
    
    cost::Float64
    c0::Vector{Float64}
    c1::Vector{Float64}
    c2::Vector{Float64}
    cost_order::Int
    source_type::String

    cg::Vector{Float64}
    slack::Int
    Ag::SparseMatrixCSC{Float64,Int64}
    Abss::SparseMatrixCSC{Float64,Int64}
    Aij::SparseMatrixCSC{Float64,Int64}
    Aji::SparseMatrixCSC{Float64,Int64}
    
    model

    PowersenseData(nbus::Int,ngen::Int,nbr::Int,nbss::Int) = new(
        nbus,ngen,nbr,nbss,zeros(nbus),zeros(nbus),zeros(ngen),zeros(ngen),zeros(ngen),zeros(ngen),zeros(nbss),zeros(nbss),zeros(nbr),
        zeros(ngen),zeros(ngen),zeros(nbss),ones(nbus),zeros(nbus), zeros(nbr),zeros(nbr),zeros(nbr),zeros(nbr),zeros(nbr),zeros(nbr),zeros(nbr),zeros(nbr),
        ones(nbus),zeros(nbus),Array{Tuple{Int64,Int64},1}(),Dict(),zeros(nbus),zeros(nbus),zeros(nbus),zeros(nbus),        
        zeros(nbus),zeros(nbus),zeros(nbr),zeros(nbr),zeros(2*nbr),zeros(2*nbr),zeros(ngen),Vector{Vector{Float64}}(),Vector{Vector{Float64}}(),
        0.0,zeros(ngen),zeros(ngen),zeros(ngen),0, " ",
        zeros(ngen),1,spzeros(nbus, ngen),spzeros(nbus, nbss),spzeros(nbus, nbr),spzeros(nbus, nbr)
    )
end

function create_PowersenseData(network::Dict{String,Any})
    m = PowersenseData(length(network["bus"]),length(network["gen"]),length(network["branch"]),length(network["sshunt"])-1);
    cost_segments = 30;
    "Assigning bss parameters"
    @simd for i=1:m.nbss sh=network["sshunt"][string(i)]; st=sh["status"]; bus=sh["shunt_bus_index"];
        @inbounds m.bmax[i] = st * sh["shmax"]; 
        @inbounds m.bmin[i] = st * sh["shmin"]; 
        m.Abss[bus,i] = 1.0;
    end
    m.slack = network["slack"];
    m.source_type=network["source_type"];
    
    m.cost_order = (m.source_type == "matpower") ? network["gen_cost_order"] : 0;
    "Assigning generator parameters"
    #@show m.ngen
    @simd for i=1:m.ngen gen=network["gen"][string(i)]; st=gen["gen_status"]; bus=gen["gen_bus_index"];
    	m.Ag[bus,i] = 1.0; 
        m.Pmin[i] = st * gen["pmin"]; 
        m.Pmax[i] = st * gen["pmax"];
        m.Qmin[i] = st * gen["qmin"]; 
        m.Qmax[i] = st * gen["qmax"]; 

        if network["source_type"] == "matpower"
        	gen = merge(gen, Dict("LinCost"=>Dict()))
        	if length(network["gen"][string(i)]["cost"]) == 3
        		m.c2[i] = gen["cost"][1];
        		m.c1[i] = gen["cost"][2];
        		m.c0[i] = gen["cost"][3];
        		gen["LinPoints"] = cost_segments;
        		l_seg = (m.Pmax[i] - m.Pmin[i])/ (gen["LinPoints"] - 1)
        		for l=1:gen["LinPoints"] 
        			p = m.Pmin[i] + (l-1) * l_seg
        			gen["LinCost"] = merge(gen["LinCost"], Dict(string(l)=>[p, m.c2[i] * p^2 + m.c1[i] * p + m.c0[i]]))
        		end
        	elseif length(network["gen"][string(i)]["cost"]) == 2
        		m.c1[i] = gen["cost"][1];
        		m.c0[i] = gen["cost"][2];
        		gen["LinPoints"] = 2;
        		gen["LinCost"] = merge(gen["LinCost"], Dict(string(1)=>[m.Pmin[i],m.c0[i]]))
        		gen["LinCost"] = merge(gen["LinCost"], Dict(string(2)=>[m.Pmax[i],m.c1[i]]))
        	elseif length(network["gen"][string(i)]["cost"]) == 1
        		m.c0[i] = gen["cost"][1];
        		gen["LinPoints"] = 1;
        		gen["LinCost"] = merge(gen["LinCost"], Dict(string(1)=>[m.Pmin[i],m.c0[i]]))
        	elseif length(network["gen"][string(i)]["cost"]) < 1
        		gen["LinPoints"] = 1;
        		gen["LinCost"] = merge(gen["LinCost"], Dict(string(1)=>[m.Pmin[i],m.c0[i]]))
        	end
        end
        
        push!(m.cgh, Float64[]); push!(m.pgh, Float64[]);
        #@show i
        m.lps[i] = gen["LinPoints"]; 
        genCost = gen["LinCost"];
        if (m.Pmin[i] <= m.Pmax[i] && !(m.Pmax[i] == 0 && m.Pmin[i] == 0))
            for l=1:m.lps[i] 
                if genCost[string(l)][2] < 10e-6
                    genCost[string(l)][2] = 10e-6
                end
                push!(m.cgh[i],genCost[string(l)][2]); 
                push!(m.pgh[i],genCost[string(l)][1]); 
            end
            if (m.Pmin[i] < m.pgh[i][1]) 
                m.pgh[i][1] = m.Pmin[i]; 
            end
            if (m.Pmax[i] > m.pgh[i][m.lps[i]]) 
                m.pgh[i][m.lps[i]] = m.Pmax[i]; 
            end
        elseif (m.Pmax[i] == 0 && m.Pmin[i] == 0) 
            m.lps[i] = 1; 
            push!(m.cgh[i],0); 
            push!(m.pgh[i],0);
        end
    end
    m.source_type = "arpae"
    "Assigning bus parameters"
    @simd for i=1:m.nbus bus=network["bus"][string(i)]; 
        m.Vmax[i] = bus["vmax"];          
        m.Vmin[i] = bus["vmin"];
        m.Pd[i] = bus["Pd"];
        m.Qd[i] = bus["Qd"];
        m.Gl[i] = bus["GL"];
        m.Bl[i] = bus["BL"];
    end
    "Assigning branch parameters"
    for i=1:m.nbr br=network["branch"][string(i)]; st=br["br_status"]; f = br["f_bus_index"]; t = br["t_bus_index"];
        g = real(1 / (br["br_r"] + br["br_x"]im));
        b = imag(1 / (br["br_r"] + br["br_x"]im))
        τ = 1/br["tap"] * cos(-br["shift"]) + (1/br["tap"] * sin(-br["shift"]))im
        α = real(τ);
        β = imag(τ);
        m.Gf[i] = st * (-g*α - b*β);
        m.Bf[i] = st * (g*β - b*α);
        m.Gt[i] = st * (-g*α + b*β);
        m.Bt[i] = st * (-g*β - b*α);
        m.gf[i] = st * (br["gfM"] + (g + br["g_fr"]) / (br["tap"] * br["tap"]));
        m.bf[i] = st * (br["bfM"] + (b + br["b_fr"]) / (br["tap"] * br["tap"]));
        m.gt[i] = st * (g + br["g_fr"]);
        m.bt[i] = st * (b + br["b_fr"]);
        m.Imax[i] = st * br["rate_a"];
        m.Aij[f,i] = 1.0;
        m.Aji[t,i] = 1.0; 
        push!(m.br, (f,t));
    end
    m.node = network["nbal"];
    return m
end



function display_preprocess_info(m::PowersenseData)
	println("---> System Information");
	println(" Buses#: ", m.nbus); 
	println(" Generators#: ", m.ngen); 
	println(" Branches#: ", m.nbr); 
	println(" Controllable Shunts#: ", m.nbss); 
	println(" Threads: ", Threads.nthreads()); 
	println(" ... time_stamp: ", time()-m.start_time, " seconds"); 
	println("\n");
end


