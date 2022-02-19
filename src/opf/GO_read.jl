function read_files_basecase(PAROTS_flag, raw_path, rop_path, inl_path, con_path)
    if (PAROTS_flag) println("Reading input DATA: "); end
    contingency=zeros(Int32,0);
    (network, bus_shift)=process_raw(PAROTS_flag, raw_path);
    network=process_inl(PAROTS_flag, inl_path, network);
    network=process_rop(PAROTS_flag, rop_path, network, bus_shift);
    #contingency=process_con(PAROTS_flag, con_path, network, bus_shift);
    #network=add_LODFs(network, contingency)
    return(network, contingency);
end

function read_files_contingency(PAROTS_flag, raw_path, rop_path, inl_path, con_path)
    if (PAROTS_flag) println("Reading input DATA: "); end
    (network, bus_shift)=process_raw(PAROTS_flag, raw_path);
    network=process_inl(PAROTS_flag, inl_path, network);
    #network=process_rop(rop_path, network, bus_shift);
    contingency=process_con(PAROTS_flag, con_path, network, bus_shift);
    return(network, contingency);
end


"--------------------------------------------------
      Readiang 'raw' file (Network Data File)
---------------------------------------------------"
function add_buses(net_bus, RawString, start, size)
    bus_index=Dict{String,Any}();
    for I=1:size
        bus=split(RawString[I+start-1], ",");
        i=parse(Int64, split(bus[1]," "; limit=0, keepempty=false)[1]);
        net_bus=merge(net_bus, Dict(string(i)=>Dict{String,Any}()));
        VM=parse(Float64, split(bus[8]," "; limit=0, keepempty=false)[1]);
        VA=parse(Float64, split(bus[9]," "; limit=0, keepempty=false)[1]);
        NVHI=parse(Float64, split(bus[10]," "; limit=0, keepempty=false)[1]);
        NVLO=parse(Float64, split(bus[11]," "; limit=0, keepempty=false)[1]);
        EVHI=parse(Float64, split(bus[12]," "; limit=0, keepempty=false)[1]);
        EVLO=parse(Float64, split(bus[13]," "; limit=0, keepempty=false)[1]);
        bus_index=merge(bus_index, Dict(string(i)=>I));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("index"=>I));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("gen"=>false));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("vm"=>VM));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("va"=>deg2rad(VA)));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("vmax"=>NVHI));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("vmin"=>NVLO));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("evmax"=>EVHI));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("evmin"=>EVLO));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("bus_type"=>1));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("Pd"=>0.0));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("Qd"=>0.0));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("GL"=>0.0));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("BL"=>0.0));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("sh_max"=>0.0));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("sh_min"=>0.0));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("bcs"=>0.0));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("gfM"=>0.0));
        net_bus[string(i)]=merge(net_bus[string(i)], Dict("bfM"=>0.0));
    end
    return(net_bus, bus_index)
end
function add_loads(net_bus, bus_index, BaseMVA, RawString, start, size)
    for i=1:size
        load=split(RawString[i+start-1], ",");
        I=parse(Int64, split(load[1]," "; limit=0, keepempty=false)[1]);
        # I=bus_index[string(I_index)];
        ST=parse(Int64, split(load[3]," "; limit=0, keepempty=false)[1]);
        PL=ST*parse(Float64, split(load[6]," "; limit=0, keepempty=false)[1])/BaseMVA+net_bus[string(I)]["Pd"];
        QL=ST*parse(Float64, split(load[7]," "; limit=0, keepempty=false)[1])/BaseMVA+net_bus[string(I)]["Qd"];
        net_bus[string(I)]=merge(net_bus[string(I)], Dict("Pd"=>PL));
        net_bus[string(I)]=merge(net_bus[string(I)], Dict("Qd"=>QL));
    end
    return(net_bus)
end

function add_shunts(net_bus, bus_index, BaseMVA, RawString, start, size)
    for i=1:size
        shunt=split(RawString[i+start-1], ",");
        I_index=parse(Int64, split(shunt[1]," "; limit=0, keepempty=false)[1]);
        I=bus_index[string(I_index)];
        ST=parse(Int64, split(shunt[3]," "; limit=0, keepempty=false)[1]);
        GL=ST*parse(Float64, split(shunt[4]," "; limit=0, keepempty=false)[1])/BaseMVA+net_bus[string(I)]["GL"];
        BL=ST*parse(Float64, split(shunt[5]," "; limit=0, keepempty=false)[1])/BaseMVA+net_bus[string(I)]["BL"];
        net_bus[string(I)]=merge(net_bus[string(I)], Dict("GL"=>GL));
        net_bus[string(I)]=merge(net_bus[string(I)], Dict("BL"=>BL));
    end
    return(net_bus)
end
function add_gens(net_bus, net_gen, bus_index, BaseMVA, RawString, start, size)
    temp_pamx_id=1; gen_buses=zeros(Int32,0);
    for i=1:size
        net_gen=merge(net_gen, Dict(string(i)=>Dict{String,Any}()));
        gen=split(RawString[i+start-1], ",");
        I=parse(Int64, split(gen[1]," "; limit=0, keepempty=false)[1]);
        I_index=bus_index[string(I)];
        ID=split(gen[2],r"\W"; limit=0, keepempty=false);
        if (length(ID)>0) ID=split(ID[1]," "; limit=0, keepempty=false)[1]; end
        St=parse(Int64, split(gen[15]," "; limit=0, keepempty=false)[1]);
        PG=parse(Float64, split(gen[3]," "; limit=0, keepempty=false)[1])/BaseMVA;
        QG=parse(Float64, split(gen[4]," "; limit=0, keepempty=false)[1])/BaseMVA;
        QT=parse(Float64, split(gen[5]," "; limit=0, keepempty=false)[1])/BaseMVA;
        QB=parse(Float64, split(gen[6]," "; limit=0, keepempty=false)[1])/BaseMVA;
        PT=parse(Float64, split(gen[17]," "; limit=0, keepempty=false)[1])/BaseMVA;
        PB=parse(Float64, split(gen[18]," "; limit=0, keepempty=false)[1])/BaseMVA;

        net_gen[string(i)]=merge(net_gen[string(i)], Dict("gen_bus_index"=>I_index));
        net_gen[string(i)]=merge(net_gen[string(i)], Dict("gen_bus"=>I));
        net_gen[string(i)]=merge(net_gen[string(i)], Dict("id"=>ID));
        net_gen[string(i)]=merge(net_gen[string(i)], Dict("gen_status"=>St));
        net_gen[string(i)]=merge(net_gen[string(i)], Dict("pg"=>PG));
        net_gen[string(i)]=merge(net_gen[string(i)], Dict("qg"=>QG));
        net_gen[string(i)]=merge(net_gen[string(i)], Dict("pmax"=>PT));
        net_gen[string(i)]=merge(net_gen[string(i)], Dict("qmax"=>QT));
        net_gen[string(i)]=merge(net_gen[string(i)], Dict("pmin"=>PB));
        net_gen[string(i)]=merge(net_gen[string(i)], Dict("qmin"=>QB));
        net_gen[string(i)]=merge(net_gen[string(i)], Dict("index"=>i));
        net_gen[string(i)]=merge(net_gen[string(i)], Dict("gov_droop"=>0.0));
        net_gen[string(i)]=merge(net_gen[string(i)], Dict("total_apf"=>0.0));

        if (St==1) net_bus[string(I)]["bus_type"]=2; net_bus[string(I)]["gen"]=true; append!(gen_buses,I); end
        if (St*PT>net_gen[string(temp_pamx_id)]["pmax"]) temp_pamx_id=i; end
    end
    net_bus[string(net_gen[string(temp_pamx_id)]["gen_bus"])]["bus_type"]=3;
    return(net_bus, net_gen, net_gen[string(temp_pamx_id)]["gen_bus"], gen_buses);
end

function add_Lbranches(net_br, bus_index, BaseMVA, RawString, start, size)
    i=0;
    for ind=1:size
        Lbr=split(RawString[ind+start-1], ",");
        St=parse(Int64, split(Lbr[14]," "; limit=0, keepempty=false)[1]);
        if (St==1) i+=1;
            net_br=merge(net_br, Dict(string(i)=>Dict{String,Any}()));
            # println("string_i: ", string(i))
            I_index=parse(Int64, split(Lbr[1]," "; limit=0, keepempty=false)[1]);
            I=bus_index[string(I_index)];
            J_index=parse(Int64, split(Lbr[2]," "; limit=0, keepempty=false)[1]);
            J=bus_index[string(J_index)];
            CKT=split(Lbr[3],r"\W"; limit=0, keepempty=false);
            if (length(CKT)>0) CKT=split(CKT[1]," "; limit=0, keepempty=false)[1]; end
            R=parse(Float64, split(Lbr[4]," "; limit=0, keepempty=false)[1]);
            X=parse(Float64, split(Lbr[5]," "; limit=0, keepempty=false)[1]);
            B=parse(Float64, split(Lbr[6]," "; limit=0, keepempty=false)[1]);
            RATEA=parse(Float64, split(Lbr[7]," "; limit=0, keepempty=false)[1])/BaseMVA;
            RATEC=parse(Float64, split(Lbr[9]," "; limit=0, keepempty=false)[1])/BaseMVA;

            net_br[string(i)]=merge(net_br[string(i)], Dict("br_r"=>R));
            net_br[string(i)]=merge(net_br[string(i)], Dict("ckt"=>CKT));
            net_br[string(i)]=merge(net_br[string(i)], Dict("br_status"=>St));
            net_br[string(i)]=merge(net_br[string(i)], Dict("br_x"=>X));
            net_br[string(i)]=merge(net_br[string(i)], Dict("b_fr"=>B/2));
            net_br[string(i)]=merge(net_br[string(i)], Dict("b_to"=>B/2));
            net_br[string(i)]=merge(net_br[string(i)], Dict("rate_a"=>RATEA));
            net_br[string(i)]=merge(net_br[string(i)], Dict("rate_c"=>RATEC));
            net_br[string(i)]=merge(net_br[string(i)], Dict("f_bus"=>I_index));
            net_br[string(i)]=merge(net_br[string(i)], Dict("t_bus"=>J_index));
            net_br[string(i)]=merge(net_br[string(i)], Dict("f_bus_index"=>I));
            net_br[string(i)]=merge(net_br[string(i)], Dict("t_bus_index"=>J));
            net_br[string(i)]=merge(net_br[string(i)], Dict("transformer"=>false));
            net_br[string(i)]=merge(net_br[string(i)], Dict("tap"=>1.0));
            net_br[string(i)]=merge(net_br[string(i)], Dict("shift"=>0.0));
            net_br[string(i)]=merge(net_br[string(i)], Dict("index"=>i));
            net_br[string(i)]=merge(net_br[string(i)], Dict("gfM"=>0.0));
            net_br[string(i)]=merge(net_br[string(i)], Dict("bfM"=>0.0));
            net_br[string(i)]=merge(net_br[string(i)], Dict("g_fr"=>0.0));
            net_br[string(i)]=merge(net_br[string(i)], Dict("g_to"=>0.0));
        end
    end
    return(net_br)
end
function add_Tbranches(net_bus, net_br, bus_index, BaseMVA, RawString, pre, start, size)
    # println("pre: ", pre);
    ind=0
    for i in 1:4:size
        Tran1=split(RawString[i+start-1], ",");
        Tran2=split(RawString[i+start], ",");
        Tran3=split(RawString[i+start+1], ",");
        Tran4=split(RawString[i+start+2], ",");
        St=parse(Int64, split(Tran1[12]," "; limit=0, keepempty=false)[1]);
        if (St==1)  ind+=1;
            net_br=merge(net_br, Dict(string(ind+pre)=>Dict{String,Any}()));
            I=parse(Int64, split(Tran1[1]," "; limit=0, keepempty=false)[1]);
            I_index=bus_index[string(I)];
            J=parse(Int64, split(Tran1[2]," "; limit=0, keepempty=false)[1]);
            J_index=bus_index[string(J)];
            CKT=split(Tran1[4],r"\W"; limit=0, keepempty=false);
            if (length(CKT)>0) CKT=split(CKT[1]," "; limit=0, keepempty=false)[1]; end
            MAG1=parse(Float64, split(Tran1[8]," "; limit=0, keepempty=false)[1]);
            MAG2=parse(Float64, split(Tran1[9]," "; limit=0, keepempty=false)[1]);

            R12=parse(Float64, split(Tran2[1]," "; limit=0, keepempty=false)[1]);
            X12=parse(Float64, split(Tran2[2]," "; limit=0, keepempty=false)[1]);
            WINDV1=parse(Float64, split(Tran3[1]," "; limit=0, keepempty=false)[1]);
            ANG1=parse(Float64, split(Tran3[3]," "; limit=0, keepempty=false)[1]);
            RATEA1=parse(Float64, split(Tran3[4]," "; limit=0, keepempty=false)[1])/BaseMVA;
            RATEC1=parse(Float64, split(Tran3[6]," "; limit=0, keepempty=false)[1])/BaseMVA;
            WINDV2=parse(Float64, split(Tran4[1]," "; limit=0, keepempty=false)[1]);

            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("br_r"=>R12));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("ckt"=>CKT));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("br_status"=>St));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("br_x"=>X12));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("b_fr"=>0.0));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("b_to"=>0.0));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("rate_a"=>RATEA1));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("rate_c"=>RATEC1));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("f_bus"=>I));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("t_bus"=>J));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("f_bus_index"=>I_index));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("t_bus_index"=>J_index));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("transformer"=>true));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("tap"=>WINDV1/WINDV2));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("shift"=>ANG1*pi/180));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("index"=>ind+pre)); #index=i
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("g_fr"=>0.0));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("g_to"=>0.0));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("gfM"=>MAG1));
            net_br[string(ind+pre)]=merge(net_br[string(ind+pre)], Dict("bfM"=>MAG2));

            gfM_total=net_bus[string(I)]["gfM"]+MAG1*St;
            bfM_total=net_bus[string(I)]["bfM"]+MAG2*St;
            net_bus[string(I)]=merge(net_bus[string(I)], Dict("gfM"=>gfM_total));
            net_bus[string(I)]=merge(net_bus[string(I)], Dict("bfM"=>bfM_total));
        end
    end
    return(net_bus,net_br)
end

function add_switch_shunts(net_bus, net_sshunt, bus_index, BaseMVA, RawString, start, size)
    net_sshunt=merge(net_sshunt, Dict("bus_index"=>zeros(Int32,0)));
    for i=1:size
        net_sshunt=merge(net_sshunt, Dict(string(i)=>Dict{String,Any}()));
        shunt=split(RawString[i+start-1], ",");
        I=parse(Int64, split(shunt[1]," "; limit=0, keepempty=false)[1]);
        I_index=bus_index[string(I)];
        STATUS=parse(Int64, split(shunt[4]," "; limit=0, keepempty=false)[1]);

        N=zeros(8); B=zeros(8); BL=zeros(8); Shmin=0; Shmax=0;
        for i=1:8
            N[i]=parse(Int64, split(shunt[2*(4+i)+1]," "; limit=0, keepempty=false)[1]);
            B[i]=parse(Float64, split(shunt[2*(5+i)]," "; limit=0, keepempty=false)[1]);
            BL[i]=N[i]*B[i]/BaseMVA; if (BL[i]>0) Shmax+=BL[i]*STATUS; elseif (BL[i]<0) Shmin+=BL[i]*STATUS; end
        end

        net_sshunt[string(i)]=merge(net_sshunt[string(i)], Dict("status"=>STATUS), Dict("shunt_bus"=>I), Dict("shunt_bus_index"=>I_index));
        net_sshunt[string(i)]=merge(net_sshunt[string(i)], Dict("shmax"=>Shmax), Dict("shmin"=>Shmin));

        Shmax+=net_bus[string(I)]["sh_max"]*STATUS;
        Shmin+=net_bus[string(I)]["sh_min"]*STATUS;
        net_bus[string(I)]=merge(net_bus[string(I)], Dict("sh_max"=>Shmax), Dict("sh_min"=>Shmin));
        push!(net_sshunt["bus_index"],I);
    end
    return(net_bus, net_sshunt);
end

function process_raw(PAROTS_flag, path_raw)
    network1=Dict{String,Any}();
    raw_string = readlines(open(path_raw));
    temp=split(raw_string[1], ","); baseMVA=parse(Float64, split(temp[2]," "; limit=0, keepempty=false)[1]);
    network1=merge(network1, Dict("baseMVA"=>baseMVA), Dict("bus"=>Dict{String,Any}()), Dict("sshunt"=>Dict{String,Any}()), Dict("gen"=>Dict{String,Any}()));
    network1=merge(network1, Dict("branch"=>Dict{String,Any}()), Dict("solution"=>Dict()));

    bus_start=4; bus_end=0; bus_check=true; bus_col=13;
    load_start=4; load_end=0; load_check=true; load_col=14;
    shunt_start=4; shunt_end=0; shunt_check=true; shunt_col=5;
    gen_start=4; gen_end=0; gen_check=true; gen_col=28;
    branch_start=4; branch_end=0; branch_check=true; branch_col=24;
    Sshunt_start=4; Sshunt_end=0; Sshunt_check=true; Sshunt_col=26;
    trans_start=4; trans_end=0; trans_check=true; Tcol1=21; Tcol2=3; Tcol3=17; Tcol4=2;
    for i=4:length(raw_string)
        temp=split(raw_string[i],","); temp_pre=split(raw_string[i-1],","); temp1=0; temp2=0; temp3=0; temp_pre_sh=0; temp_sh=0;
        if (length(temp)>=2) temp_sh=split(temp[2]," "; limit=0, keepempty=false)[1]; temp_sh=split(temp[2],""; limit=0, keepempty=false)[1]; end
        if (length(temp_pre)>=2) temp_pre_sh=split(temp_pre[2]," "; limit=0, keepempty=false)[1]; temp_pre_sh=split(temp_pre[2],""; limit=0, keepempty=false)[1]; end
        if (i<length(raw_string)-3) temp1=split(raw_string[i+1],","); temp2=split(raw_string[i+2],","); temp3=split(raw_string[i+3],","); end

        if (length(temp)==bus_col && bus_check) bus_start=i; bus_check=false; end
        if (length(temp_pre)==bus_col && length(temp)!=bus_col) bus_end=i-1; end

        if (length(temp)==load_col && load_check) load_start=i; load_check=false; end
        if (length(temp_pre)==load_col && length(temp)!=load_col) load_end=i-1; end

        if (length(temp)==shunt_col && temp_sh=="'" && shunt_check) shunt_start=i; shunt_check=false; end
        if (length(temp_pre)==shunt_col && (length(temp)!=shunt_col || temp_sh!="'") && temp_pre_sh=="'") shunt_end=i-1; end

        if (length(temp)==gen_col && gen_check) gen_start=i; gen_check=false; end
        if (length(temp_pre)==gen_col && length(temp)!=gen_col) gen_end=i-1; end

        if (length(temp)==branch_col && branch_check) branch_start=i; branch_check=false; end
        if (length(temp_pre)==branch_col && length(temp)!=branch_col) branch_end=i-1; end

        if (length(temp)==Tcol1 && length(temp1)==Tcol2 && length(temp2)==Tcol3 && length(temp3)==Tcol4 && trans_check) trans_start=i; trans_check=false; end
        if (length(temp)==Tcol4 && length(temp_pre)==Tcol3 && length(temp1)!=Tcol1) trans_end=i; end

        if (length(temp)==Sshunt_col && Sshunt_check) Sshunt_start=i; Sshunt_check=false; end
        if (length(temp_pre)==Sshunt_col && length(temp)!=Sshunt_col) Sshunt_end=i-1; end
    end

    (network1["bus"], bus_shift)=add_buses(network1["bus"], raw_string, bus_start, bus_end-bus_start+1);
    network1["bus"]=add_loads(network1["bus"], bus_shift, network1["baseMVA"], raw_string, load_start, load_end-load_start+1);
    network1["bus"]=add_shunts(network1["bus"], bus_shift, network1["baseMVA"], raw_string, shunt_start, shunt_end-shunt_start+1);
    (network1["bus"], network1["gen"], slack, gen_buses) =add_gens(network1["bus"], network1["gen"], bus_shift, network1["baseMVA"], raw_string, gen_start, gen_end-gen_start+1);
    network1=merge(network1, Dict("gen_buses"=>gen_buses, "slack"=>slack));
    network1["branch"]=add_Lbranches(network1["branch"], bus_shift, network1["baseMVA"], raw_string, branch_start, branch_end-branch_start+1);
    (network1["bus"], network1["branch"])=add_Tbranches(network1["bus"], network1["branch"], bus_shift, network1["baseMVA"], raw_string, length(network1["branch"]), trans_start, trans_end-trans_start+1);
    (network1["bus"], network1["sshunt"])=add_switch_shunts(network1["bus"], network1["sshunt"], bus_shift, network1["baseMVA"], raw_string, Sshunt_start, Sshunt_end-Sshunt_start+1);
    network1=merge(network1, Dict("nbal"=>Dict{String,Any}(), "source_type"=>"arpae"));
    for i=1:length(network1["bus"])
        network1["nbal"]=merge(network1["nbal"], Dict(string(i)=>Dict()));
        network1["nbal"][string(i)]=merge(network1["nbal"][string(i)], Dict("gen"=>zeros(Int32,0)), Dict("Sij"=>zeros(Int32,0)), Dict("Sji"=>zeros(Int32,0)), Dict("shunt"=>zeros(Int32,0)));
    end

    for g=1:length(network1["gen"]) gen=network1["gen"][string(g)];
        if (gen["gen_bus_index"]*gen["gen_status"]>0) push!(network1["nbal"][string(gen["gen_bus_index"])]["gen"],g); end
    end
    # println("nbr: ", length(network1["branch"]));
    for b=1:length(network1["branch"]) br=network1["branch"][string(b)];
        if (br["f_bus_index"]*br["br_status"]>0) push!(network1["nbal"][string(br["f_bus_index"])]["Sij"], b); end
        if (br["t_bus_index"]*br["br_status"]>0) push!(network1["nbal"][string(br["t_bus_index"])]["Sji"], b); end
    end
    for sh=1:length(network1["sshunt"])-1 shunt=network1["sshunt"][string(sh)];
        if (shunt["shunt_bus_index"]*shunt["status"]>0) push!(network1["nbal"][string(shunt["shunt_bus_index"])]["shunt"],sh); end
    end
    # if (PAROTS_flag>0) println("The RAW file was read successfully ...") end
    return(network1, bus_shift)
end

"--------------------------------------------------
   Readiang 'rop' file (Generator Cost Data File)
---------------------------------------------------"
function Generator_Dispatch_data(rop_string, gen_index, bus_index, GDD_in_start, size, generator)
    for i=1:size
        temp = split(rop_string[GDD_in_start+i-1],",");
        temp[2]=split(temp[2],r"\W"; limit=0, keepempty=false)[1];
        ID=split(temp[2]," "; limit=0, keepempty=false)[1];
        bus_no=parse(Int64,temp[1]);
        bus_no_index=bus_index[string(bus_no)];
        DSPTBL=parse(Int64,temp[4]);
        generator[string(gen_index[string(DSPTBL)])]=merge(generator[string(gen_index[string(DSPTBL)])], Dict("GenID"=>ID), Dict("bus"=>bus_no), Dict("bus_index"=>bus_no_index));
    end
    return(generator)
end

function process_rop(PAROTS_flag, path_rop, network, bus_index)
    rop_string = readlines(open(path_rop));
    generator=Dict{String,Any}();
    "Generator Dispatch Data"
    GDU_start=0; GDU_check=true; GDU_col=4; GDU_end=0;
    APDT_start=0; APDT_check=true; APDT_col=7; #APDT_end=0;
    PLCT_start=0; PLCT_check=true; PLCT_col=3; #PLCT_end=0;

    for i=1:length(rop_string)
        temp=split(rop_string[i],",");
        temp_pre= (i>1)  ? split(rop_string[i-1],",") : ["0 "," 0"] ;
	       temp_1=split(temp[1]," "; limit=0, keepempty=false)[1];
	          temp_1_pre=split(temp_pre[1]," "; limit=0, keepempty=false)[1];
        if (length(temp)==GDU_col && GDU_check && temp_1!="0") GDU_start=i; GDU_check=false; end
        if (length(temp_pre)==GDU_col && length(temp)!=GDU_col && temp_1_pre!="0") GDU_end=i-1; end
        if (length(temp)==APDT_col && APDT_check && temp_1!="0") APDT_start=i; APDT_check=false; end
        if (length(temp)==PLCT_col && PLCT_check && temp_1!="0") PLCT_start=i; PLCT_check=false; end
    end
    size=GDU_end-GDU_start+1;
    NPAIRSit=PLCT_start;
    NPAIRS=0;
    gen_index=Dict{String,Any}();
    for i=1:size
        temp = split(rop_string[NPAIRSit],","); NPAIRS = parse(Int64, temp[3]); nTable=parse(Int64, temp[1]);
        temp=split(rop_string[NPAIRSit+1],","); X=parse(Float64, temp[1]); Y=parse(Float64, temp[2]);
        generator=merge(generator, Dict(string(nTable)=>Dict("LinPoints"=>NPAIRS, "LinCost"=>Dict())));
        baseMVA=network["baseMVA"];
        generator[string(nTable)]["LinCost"]=merge(generator[string(nTable)]["LinCost"], Dict("1"=>[X/baseMVA,Y/baseMVA]));
        for c=2:NPAIRS
            temp=split(rop_string[NPAIRSit+c],","); X=parse(Float64, temp[1]); Y=parse(Float64, temp[2]);
            generator[string(nTable)]["LinCost"]=merge(generator[string(nTable)]["LinCost"], Dict(string(c)=>[X/baseMVA,Y/baseMVA]));
        end
        NPAIRSit+=NPAIRS+1;
    end
    for i=1:size
        temp1 = split(rop_string[APDT_start+i-1],",");
        index = parse(Int64, temp1[7]); TBL = parse(Int64, temp1[1]); gen_index=merge(gen_index, Dict(string(TBL)=>index));
        generator[string(index)]=merge(generator[string(index)], Dict("gov_droop"=>0.0));
    end
    for i=1:size
        temp2 = split(rop_string[GDU_start+i-1],",");
        temp2[2]=split(temp2[2],r"\W"; limit=0, keepempty=false)[1];
        ID=split(temp2[2]," "; limit=0, keepempty=false)[1];
        bus_no=parse(Int64,temp2[1]);
        bus_no_index=bus_index[string(bus_no)];
        DSPTBL=parse(Int64,temp2[4]);
        generator[string(gen_index[string(DSPTBL)])]=merge(generator[string(gen_index[string(DSPTBL)])], Dict("GenID"=>ID), Dict("bus"=>bus_no), Dict("bus_index"=>bus_no_index));
    end
    for i=1:length(generator) if (!haskey(generator[string(i)], "LinPoints") && (PAROTS_flag)) println("generator missing info at ", i) end end
    if (size==length(network["gen"])) Gen_Dict=deepcopy(network["gen"]);
        for i=1:size
            for j in Gen_Dict
                bg=j[2]["gen_bus"]; Ug=j[2]["id"];
                if (generator[string(i)]["bus"]==bg && generator[string(i)]["GenID"]==Ug)
                    LC=generator[string(i)]["LinCost"]; LP=generator[string(i)]["LinPoints"];
                    network["gen"][string(j[1])]=merge(network["gen"][string(j[1])], Dict("LinCost"=>LC), Dict("LinPoints"=>LP));
                    delete!(Gen_Dict, j[1]);
                end
            end
        end
        # if (PAROTS_flag>0) println("The ROP file was read successfully ...") end
    end
    # elseif (PAROTS_flag) println("ERROR | The number of generators in RAW and ROP files don't match") end
    return(network)
end

"--------------------------------------------------
Readiang 'inl' file (Governor Response Data File)
---------------------------------------------------"
function process_inl(PAROTS_flag, path_inl, network)
    inl_string = readlines(open(path_inl));
    Gen_Dict=deepcopy(network["gen"]);
    for i=1:length(inl_string)
        temp=split(inl_string[i],",");
        if (length(temp)==7)
            bus=parse(Int64, split(temp[1]," "; limit=0, keepempty=false)[1]);
            ID_temp=split(temp[2],r"\W"; limit=0, keepempty=false)[1]; genID=split(ID_temp," "; limit=0, keepempty=false)[1];
            for j in Gen_Dict
                bg=j[2]["gen_bus"]; Ug=j[2]["id"];
                if (bus==bg && genID==Ug)
                    gov_droop=parse(Float64, temp[6]);
                    network["gen"][string(j[1])]=merge(network["gen"][string(j[1])], Dict("gov_droop"=>gov_droop));
                    delete!(Gen_Dict, j[1]);
                end
            end
        end
    end
    sum_apf=0
    for i=1:length(network["gen"]) gen=network["gen"][string(i)]; st=gen["gen_status"]
        sum_apf+=st*gen["gov_droop"];
    end
    for i=1:length(network["gen"]) gen=network["gen"][string(i)]; st=gen["gen_status"]
        total_apf=sum_apf-st*gen["gov_droop"];
        network["gen"][string(i)]=merge(network["gen"][string(i)], Dict("total_apf"=>total_apf));
    end
    # if (PAROTS_flag>0) println("The INL file was read successfully ..."); end
    return(network);
end

"-------------------------------------------------------
Readiang 'con' file (Contingency Description Data File)
--------------------------------------------------------"
function adding_con_line_info(contingency, counts, Branch_Dict)
    fc=contingency[string(counts)]["f_bus"];
    tc=contingency[string(counts)]["t_bus"];
    CTc=contingency[string(counts)]["CT"];
    for j in Branch_Dict
        fb=j[2]["f_bus"];
        tb=j[2]["t_bus"];
        CTb=j[2]["ckt"];
        if (fb==fc && tb==tc && CTb==CTc)
            index=j[2]["index"];
            contingency[string(counts)]=merge(contingency[string(counts)], Dict("branch"=>index));
        end
    end
    return(contingency)
end

function adding_con_gen_info(contingency, counts, contingency_G, Gcounts, Gen_Dict)
    bc=contingency[string(counts)]["bus"]; Uc=contingency[string(counts)]["unit"];
    for j in Gen_Dict
        bg=j[2]["gen_bus"]; Ug=j[2]["id"];
        if (bc==bg && Uc==Ug)
            contingency[string(counts)]=merge(contingency[string(counts)], Dict("GenNum"=>j[2]["index"]));
            contingency_G[string(Gcounts)]=merge(contingency_G[string(Gcounts)], Dict("GenNum"=>j[2]["index"]));
        end
    end
    return(contingency, contingency_G);
end

function process_con(PAROTS_flag, path_con, network, bus_index)
    con_string = readlines(open(path_con));
    contingency=Dict{String,Any}(); contingency_G=Dict{String,Any}();
    counts=0; Gcounts=0
    for i=1:length(con_string)-1
        temp1 = split(con_string[i]," "; limit=0, keepempty=false); temp2 = split(con_string[i+1]," "; limit=0, keepempty=false);

        if (length(temp1)==2 && length(temp2)<6)
            label=temp1[2]; counts+=1;
            contingency[string(counts)]=Dict("status"=>false);
            contingency[string(counts)]=merge(contingency[string(counts)], Dict("label"=>label), Dict("type"=>NaN));
        elseif (length(temp1)==2 && length(temp2)==10)
            label=temp1[2]; type = split(temp1[2],""; limit=0, keepempty=false)[1]; counts+=1; f_bus_index=parse(Int64,temp2[5]); f_bus=bus_index[string(f_bus_index)];
            contingency[string(counts)]=Dict("type"=>type, "status"=>true, "label"=>label, "f_bus"=>f_bus_index, "f_bus_index"=>f_bus);
            t_bus_index=parse(Int64,temp2[8]); t_bus=bus_index[string(t_bus_index)]; CT_bus=temp2[10];
            contingency[string(counts)]=merge(contingency[string(counts)], Dict("t_bus"=>t_bus_index, "t_bus_index"=>t_bus, "CT"=>CT_bus));
            contingency=adding_con_line_info(contingency, counts, network["branch"]);
        elseif (length(temp1)==2 && length(temp2)==6)
            label=temp1[2];
            counts+=1;
            Gcounts+=1;
            unit=temp2[3];
            f_bus_index=parse(Int64,temp2[6]);
            f_bus=bus_index[string(f_bus_index)];
            contingency[string(counts)]=Dict("type"=>"G", "status"=>true, "label"=>label, "unit"=>unit, "bus"=>f_bus_index, "bus_index"=>f_bus);
            contingency_G[string(Gcounts)]=Dict("type"=>"G", "label"=>label, "unit"=>unit, "bus"=>f_bus_index, "bus_index"=>f_bus);
            (contingency, contingency_G)=adding_con_gen_info(contingency, counts, contingency_G, Gcounts, network["gen"]);
        end
    end

    # if (PAROTS_flag>0) println("The CON file was read successfully ..."); end
    return(contingency);
end

"--------------------------------------------------
   Readiang 'raw' file for LODF calculation
---------------------------------------------------"
function process_net_LODF(PAROTS_flag, path_raw)
    network1=Dict{String,Any}();
    raw_string = readlines(open(path_raw));
    temp=split(raw_string[1], ","); BaseMVA=parse(Float64, split(temp[2]," "; limit=0, keepempty=false)[1]);
    network1=merge(network1, Dict("baseMVA"=>BaseMVA), Dict("bus"=>Dict{String,Any}()), Dict("sshunt"=>Dict{String,Any}()), Dict("gen"=>Dict{String,Any}()));
    network1=merge(network1, Dict("branch"=>Dict{String,Any}()), Dict("solution"=>Dict()));

    bus_start=4; bus_end=0; bus_check=true; bus_col=13;
    gen_start=4; gen_end=0; gen_check=true; gen_col=28;
    branch_start=4; branch_end=0; branch_check=true; branch_col=24;
    trans_start=4; trans_end=0; trans_check=true; Tcol1=21; Tcol2=3; Tcol3=17; Tcol4=2;
    for i=4:length(raw_string)
        temp=split(raw_string[i],","); temp_pre=split(raw_string[i-1],","); temp1=0; temp2=0; temp3=0; temp_pre_sh=0; temp_sh=0;
        if (length(temp)>=2) temp_sh=split(temp[2]," "; limit=0, keepempty=false)[1]; temp_sh=split(temp[2],""; limit=0, keepempty=false)[1]; end
        if (length(temp_pre)>=2) temp_pre_sh=split(temp_pre[2]," "; limit=0, keepempty=false)[1]; temp_pre_sh=split(temp_pre[2],""; limit=0, keepempty=false)[1]; end
        if (i<length(raw_string)-3) temp1=split(raw_string[i+1],","); temp2=split(raw_string[i+2],","); temp3=split(raw_string[i+3],","); end

        if (length(temp)==bus_col && bus_check) bus_start=i; bus_check=false; end
        if (length(temp_pre)==bus_col && length(temp)!=bus_col) bus_end=i-1; end

        if (length(temp)==gen_col && gen_check) gen_start=i; gen_check=false; end
        if (length(temp_pre)==gen_col && length(temp)!=gen_col) gen_end=i-1; end

        if (length(temp)==branch_col && branch_check) branch_start=i; branch_check=false; end
        if (length(temp_pre)==branch_col && length(temp)!=branch_col) branch_end=i-1; end

        if (length(temp)==Tcol1 && length(temp1)==Tcol2 && length(temp2)==Tcol3 && length(temp3)==Tcol4 && trans_check) trans_start=i; trans_check=false; end
        if (length(temp)==Tcol4 && length(temp_pre)==Tcol3 && length(temp1)!=Tcol1) trans_end=i; end
    end
    bus_size=bus_end-bus_start+1; bus_shift=Dict{String,Any}();
    for i=1:bus_size
        bus=split(raw_string[i+bus_start-1], ",");
        I=parse(Int64, split(bus[1]," "; limit=0, keepempty=false)[1]);
        bus_shift=merge(bus_shift, Dict(string(I)=>i));
    end
    gen_size=gen_end-gen_start+1;
    temp_pamx_id=1; max_gen=0.0;
    for i=1:gen_size
        gen=split(raw_string[i+gen_start-1], ",");
        I_index=parse(Int64, split(gen[1]," "; limit=0, keepempty=false)[1]);
        I=bus_shift[string(I_index)];
        St=parse(Int64, split(gen[15]," "; limit=0, keepempty=false)[1]);
        PT=parse(Float64, split(gen[17]," "; limit=0, keepempty=false)[1])/BaseMVA;
        if (max_gen<St*PT) max_gen=St*PT; temp_pamx_id=i; end
    end
    network1=merge(network1, Dict("slack"=>temp_pamx_id));
    network1=merge(network1, Dict("nbus"=>bus_size));
    branch_size=branch_end-branch_start+1; i=0;
    for ind=1:branch_size
        Lbr=split(raw_string[ind+branch_start-1], ",");
        St=parse(Int64, split(Lbr[14]," "; limit=0, keepempty=false)[1]);
        if (St==1) i+=1;
            network1["branch"]=merge(network1["branch"], Dict(string(i)=>Dict{String,Any}()));
            I_index=parse(Int64, split(Lbr[1]," "; limit=0, keepempty=false)[1]);
            I=bus_shift[string(I_index)];
            J_index=parse(Int64, split(Lbr[2]," "; limit=0, keepempty=false)[1]);
            J=bus_shift[string(J_index)];
            CKT=split(Lbr[3],r"\W"; limit=0, keepempty=false);
            if (length(CKT)>0) CKT=split(CKT[1]," "; limit=0, keepempty=false)[1]; end

            R=parse(Float64, split(Lbr[4]," "; limit=0, keepempty=false)[1]);
            X=parse(Float64, split(Lbr[5]," "; limit=0, keepempty=false)[1]);
            B=parse(Float64, split(Lbr[6]," "; limit=0, keepempty=false)[1]);
            RATEC=parse(Float64, split(Lbr[9]," "; limit=0, keepempty=false)[1])/BaseMVA;

            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("br_r"=>R));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("br_status"=>St, "index"=>i));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("br_x"=>X,"ratec"=>RATEC));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("b_fr"=>B/2));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("b_to"=>B/2));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("f_bus"=>I_index, "t_bus"=>J_index, "ckt"=>CKT));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("f_bus_index"=>I));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("t_bus_index"=>J));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("tap"=>1.0));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("shift"=>0.0));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("gfM"=>0.0));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("bfM"=>0.0));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("g_fr"=>0.0));
            network1["branch"][string(i)]=merge(network1["branch"][string(i)], Dict("g_to"=>0.0));
        end
    end
    branch_size=length(network1["branch"]);
    tbranch_size=trans_end-trans_start+1; ind=0;
    for i in 1:4:tbranch_size
        Tran1=split(raw_string[i+trans_start-1], ",");
        Tran2=split(raw_string[i+trans_start], ",");
        Tran3=split(raw_string[i+trans_start+1], ",");
        Tran4=split(raw_string[i+trans_start+2], ",");
        St=parse(Int64, split(Tran1[12]," "; limit=0, keepempty=false)[1]);
        if (St==1) ind+=1;
            network1["branch"]=merge(network1["branch"], Dict(string(ind+branch_size)=>Dict{String,Any}()));
            I_index=parse(Int64, split(Tran1[1]," "; limit=0, keepempty=false)[1]);
            I=bus_shift[string(I_index)];
            J_index=parse(Int64, split(Tran1[2]," "; limit=0, keepempty=false)[1]);
            J=bus_shift[string(J_index)];
            CKT=split(Tran1[4],r"\W"; limit=0, keepempty=false);
            if (length(CKT)>0) CKT=split(CKT[1]," "; limit=0, keepempty=false)[1]; end
            MAG1=parse(Float64, split(Tran1[8]," "; limit=0, keepempty=false)[1]);
            MAG2=parse(Float64, split(Tran1[9]," "; limit=0, keepempty=false)[1]);
            R12=parse(Float64, split(Tran2[1]," "; limit=0, keepempty=false)[1]);
            X12=parse(Float64, split(Tran2[2]," "; limit=0, keepempty=false)[1]);
            WINDV1=parse(Float64, split(Tran3[1]," "; limit=0, keepempty=false)[1]);
            ANG1=parse(Float64, split(Tran3[3]," "; limit=0, keepempty=false)[1]);
            WINDV2=parse(Float64, split(Tran4[1]," "; limit=0, keepempty=false)[1]);
            RATEC1=parse(Float64, split(Tran3[6]," "; limit=0, keepempty=false)[1])/BaseMVA;

            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("br_r"=>R12));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("br_status"=>St, "index"=>ind+branch_size));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("br_x"=>X12));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("b_fr"=>0.0,"ratec"=>RATEC1));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("b_to"=>0.0));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("f_bus"=>I_index, "t_bus"=>J_index, "ckt"=>CKT));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("f_bus_index"=>I));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("t_bus_index"=>J));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("transformer"=>true));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("tap"=>WINDV1/WINDV2));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("shift"=>deg2rad(ANG1)));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("g_fr"=>0.0));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("g_to"=>0.0));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("gfM"=>MAG1));
            network1["branch"][string(ind+branch_size)]=merge(network1["branch"][string(ind+branch_size)], Dict("bfM"=>MAG2));
        end
    end
    return(network1, bus_shift)
end
"--------------------------------------------------
        find slack bus from 'raw' file
---------------------------------------------------"
function find_slack(path_raw)
    raw_string = readlines(open(path_raw));
    temp=split(raw_string[1], ","); BaseMVA=parse(Float64, split(temp[2]," "; limit=0, keepempty=false)[1]);
    bus_start=4; bus_end=0; bus_check=true; bus_col=13;
    gen_start=4; gen_end=0; gen_check=true; gen_col=28;
    for i=4:length(raw_string)
        temp=split(raw_string[i],","); temp_pre=split(raw_string[i-1],",");
        if (length(temp)==bus_col && bus_check) bus_start=i; bus_check=false; end
        if (length(temp_pre)==bus_col && length(temp)!=bus_col) bus_end=i-1; end
        if (length(temp)==gen_col && gen_check) gen_start=i; gen_check=false; end
        if (length(temp_pre)==gen_col && length(temp)!=gen_col) gen_end=i-1; end
    end
    bus_size=bus_end-bus_start+1; bus_shift=Dict{String,Any}();
    for i=1:bus_size
        bus=split(raw_string[i+bus_start-1], ",");
        I=parse(Int64, split(bus[1]," "; limit=0, keepempty=false)[1]);
        bus_shift=merge(bus_shift, Dict(string(I)=>i));
    end
    gen_size=gen_end-gen_start+1; index=0;
    for i=1:gen_size
        gen=split(raw_string[i+gen_start-1], ",");
        I_index=parse(Int64, split(gen[1]," "; limit=0, keepempty=false)[1]);
        I=bus_shift[string(I_index)];
        St=parse(Int64, split(gen[15]," "; limit=0, keepempty=false)[1]);
        if (St==1) index=I; break; end
    end
    return(index)
end

"--------------------------------------------------
 Readiang 'raw' file for node4 (Network Data File)
---------------------------------------------------"
function process_raw_node4(PAROTS_flag, path_raw)
    network1=Dict{String,Any}();
    raw_string = readlines(open(path_raw));
    temp=split(raw_string[1], ","); baseMVA=parse(Float64, split(temp[2]," "; limit=0, keepempty=false)[1]);
    network1=merge(network1, Dict("baseMVA"=>baseMVA), Dict("bus"=>Dict{String,Any}()), Dict("sshunt"=>Dict{String,Any}()), Dict("gen"=>Dict{String,Any}()));
    network1=merge(network1, Dict("branch"=>Dict{String,Any}()), Dict("solution"=>Dict()));

    bus_start=4; bus_end=0; bus_check=true; bus_col=13;
    load_start=4; load_end=0; load_check=true; load_col=14;
    shunt_start=4; shunt_end=0; shunt_check=true; shunt_col=5;
    gen_start=4; gen_end=0; gen_check=true; gen_col=28;
    branch_start=4; branch_end=0; branch_check=true; branch_col=24;
    Sshunt_start=4; Sshunt_end=0; Sshunt_check=true; Sshunt_col=26;
    trans_start=4; trans_end=0; trans_check=true; Tcol1=21; Tcol2=3; Tcol3=17; Tcol4=2;
    for i=4:length(raw_string)
        temp=split(raw_string[i],","); temp_pre=split(raw_string[i-1],","); temp1=0; temp2=0; temp3=0; temp_pre_sh=0; temp_sh=0;
        if (length(temp)>=2) temp_sh=split(temp[2]," "; limit=0, keepempty=false)[1]; temp_sh=split(temp[2],""; limit=0, keepempty=false)[1]; end
        if (length(temp_pre)>=2) temp_pre_sh=split(temp_pre[2]," "; limit=0, keepempty=false)[1]; temp_pre_sh=split(temp_pre[2],""; limit=0, keepempty=false)[1]; end
        if (i<length(raw_string)-3) temp1=split(raw_string[i+1],","); temp2=split(raw_string[i+2],","); temp3=split(raw_string[i+3],","); end

        if (length(temp)==bus_col && bus_check) bus_start=i; bus_check=false; end
        if (length(temp_pre)==bus_col && length(temp)!=bus_col) bus_end=i-1; end

        if (length(temp)==load_col && load_check) load_start=i; load_check=false; end
        if (length(temp_pre)==load_col && length(temp)!=load_col) load_end=i-1; end

        if (length(temp)==shunt_col && temp_sh=="'" && shunt_check) shunt_start=i; shunt_check=false; end
        if (length(temp_pre)==shunt_col && (length(temp)!=shunt_col || temp_sh!="'") && temp_pre_sh=="'") shunt_end=i-1; end

        if (length(temp)==gen_col && gen_check) gen_start=i; gen_check=false; end
        if (length(temp_pre)==gen_col && length(temp)!=gen_col) gen_end=i-1; end

        if (length(temp)==branch_col && branch_check) branch_start=i; branch_check=false; end
        if (length(temp_pre)==branch_col && length(temp)!=branch_col) branch_end=i-1; end

        if (length(temp)==Tcol1 && length(temp1)==Tcol2 && length(temp2)==Tcol3 && length(temp3)==Tcol4 && trans_check) trans_start=i; trans_check=false; end
        if (length(temp)==Tcol4 && length(temp_pre)==Tcol3 && length(temp1)!=Tcol1) trans_end=i; end

        if (length(temp)==Sshunt_col && Sshunt_check) Sshunt_start=i; Sshunt_check=false; end
        if (length(temp_pre)==Sshunt_col && length(temp)!=Sshunt_col) Sshunt_end=i-1; end
    end

    (network1["bus"], bus_shift)=add_buses(network1["bus"], raw_string, bus_start, bus_end-bus_start+1);
    (network1["bus"], network1["gen"], slack) =add_gens(network1["bus"], network1["gen"], bus_shift, network1["baseMVA"], raw_string, gen_start, gen_end-gen_start+1);
    # @save "slack.jld2" slack
    network1["bus"]=add_loads(network1["bus"], bus_shift, network1["baseMVA"], raw_string, load_start, load_end-load_start+1);
    network1["bus"]=add_shunts(network1["bus"], bus_shift, network1["baseMVA"], raw_string, shunt_start, shunt_end-shunt_start+1);
    network1=merge(network1, Dict("slack"=>slack));
    network1["branch"]=add_Lbranches(network1["branch"], bus_shift, network1["baseMVA"], raw_string, branch_start, branch_end-branch_start+1);
    (network1["bus"], network1["branch"])=add_Tbranches(network1["bus"], network1["branch"], bus_shift, network1["baseMVA"], raw_string, branch_end-branch_start+1, trans_start, trans_end-trans_start+1);
    (network1["bus"], network1["sshunt"])=add_switch_shunts(network1["bus"], network1["sshunt"], bus_shift, network1["baseMVA"], raw_string, Sshunt_start, Sshunt_end-Sshunt_start+1);
    network1=merge(network1, Dict("nbal"=>Dict{String,Any}()));
    for i=1:length(network1["bus"])
        network1["nbal"]=merge(network1["nbal"], Dict(string(i)=>Dict()));
        network1["nbal"][string(i)]=merge(network1["nbal"][string(i)], Dict("gen"=>zeros(Int32,0)), Dict("Sij"=>zeros(Int32,0)), Dict("Sji"=>zeros(Int32,0)), Dict("shunt"=>zeros(Int32,0)));
    end

    for g=1:length(network1["gen"]) gen=network1["gen"][string(g)];
        if (gen["gen_bus_index"]*gen["gen_status"]>0) push!(network1["nbal"][string(gen["gen_bus_index"])]["gen"],g); end
    end
    for b=1:length(network1["branch"]) br=network1["branch"][string(b)];
        if (br["f_bus_index"]*br["br_status"]>0) push!(network1["nbal"][string(br["f_bus_index"])]["Sij"], b); end
        if (br["t_bus_index"]*br["br_status"]>0) push!(network1["nbal"][string(br["t_bus_index"])]["Sji"], b); end
    end
    for sh=1:length(network1["sshunt"])-1 shunt=network1["sshunt"][string(sh)];
        if (shunt["shunt_bus_index"]*shunt["status"]>0) push!(network1["nbal"][string(shunt["shunt_bus_index"])]["shunt"],sh); end
    end
    return(network1, bus_shift)
end

"--------------------------------------------------
 Readiang 'con' file for node5 (Network Data File)
---------------------------------------------------"
function adding_con_line_info_node5(contingency, Branch_Dict, LconList)
    fc=contingency["f_bus"];
    tc=contingency["t_bus"];
    CTc=contingency["CT"];
    for j in Branch_Dict
        fb=j[2]["f_bus"];
        tb=j[2]["t_bus"];
        CTb=j[2]["ckt"];
        if (fb==fc && tb==tc && CTb==CTc)
            index=j[2]["index"];
            append!(LconList,index);
        end
    end
    return(LconList)
end
function process_con_node5(PAROTS_flag, path_con, branch, bus_index)
    con_string = readlines(open(path_con));
    LconList=zeros(Int32,0);
    for i=1:length(con_string)-1
        temp1 = split(con_string[i]," "; limit=0, keepempty=false); temp2 = split(con_string[i+1]," "; limit=0, keepempty=false);
        if (length(temp1)==2 && length(temp2)==10)
            f_bus_index=parse(Int64,temp2[5]); f_bus=bus_index[string(f_bus_index)];
            t_bus_index=parse(Int64,temp2[8]); t_bus=bus_index[string(t_bus_index)]; CT_bus=temp2[10];
            contingency=Dict("f_bus"=>f_bus_index, "t_bus"=>t_bus_index, "CT"=>CT_bus);
            LconList=adding_con_line_info_node5(contingency, branch, LconList);
        end
    end
    # data=Dict("data"=>LconList);
    return(LconList)
    # @save "LconList.bson" data
    # @save "LconList.bson" Dict("data"=>LconList)
end
