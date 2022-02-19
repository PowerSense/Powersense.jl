include("type.jl")

function create_opf_model!(data::PowersenseData;
    formulation = PBRAPVmodel,
    initialize = true,
    FACTS_bShunt = true,
    FACTS_bSeries = true,
    box_constraints = false,
    obj_type::Symbol = :linear)    #obj_type = "linear" or "quadratic"

    "JuMP Model for ACOPF formulations"
    data.model = JuMP.Model();

    if formulation ∈ [PBDCmodel, PNDCmodel, DCmodel]
        create_dcopf_model!(data, formulation = formulation, initialize = initialize, FACTS_bSeries = FACTS_bSeries, obj_type = obj_type);
        return
    end

    #Variable initialization
    if initialize
        #Adding variables: modeling voltages for different ACOPF formulations 
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel] 
            data._θ=JuMP.@variable(data.model, θ[i=1:data.nbus], start=data.θ[i]);                       
            data._V=JuMP.@variable(data.model, V[i=1:data.nbus], start=data.V[i]);
        elseif formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel, PBRAWVmodel, CBRAWVmodel, PNRAWVmodel]
            data._vi=JuMP.@variable(data.model, vi[i=1:data.nbus], start=data.vi[i]);                       
            data._vr=JuMP.@variable(data.model, vr[i=1:data.nbus], start=data.vr[i]);
        end
        if formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            data._Wd=JuMP.@variable(data.model, Wd[i=1:data.nbus], start=data.Wd[i]);        
        end
        if formulation ∈ [PBRAWVmodel, PNRAWVmodel] 
            data._Wr=JuMP.@variable(data.model, Wr[i=1:data.nbr], start=data.Wr[i]);        
            data._Wi=JuMP.@variable(data.model, Wi[i=1:data.nbr], start=data.Wi[i]);
        end

        #Adding variables: modeling generator injections for different ACOPF formulations
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
                            PBRARVmodel, CBRARVmodel, PNRARVmodel, 
                            PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            data._pg=JuMP.@variable(data.model, pg[i=1:data.ngen], start=data.Pg[i]);        
            data._qg=JuMP.@variable(data.model, qg[i=1:data.ngen], start=data.Qg[i]);
        end

        #Adding variables: modeling branch flows for different ACOPF formulations
        if formulation ∈ [PBRAPVmodel, PBRARVmodel, PBRAWVmodel] 
            data._Pij=JuMP.@variable(data.model, Pij[i=1:data.nbr], start=data.Pij[i]);                      
            data._Qij=JuMP.@variable(data.model, Qij[i=1:data.nbr], start=data.Qij[i]);
            data._Pji=JuMP.@variable(data.model, Pji[i=1:data.nbr], start=data.Pji[i]);                     
            data._Qji=JuMP.@variable(data.model, Qji[i=1:data.nbr], start=data.Qji[i]);
        end
        if formulation ∈ [CBRARVmodel, CBRAWVmodel] 
            data._Irij=JuMP.@variable(data.model, Irij[i=1:data.nbr], start=data.Irij[i]);                      
            data._Iiij=JuMP.@variable(data.model, Iiij[i=1:data.nbr], start=data.Iiij[i]);
            data._Irji=JuMP.@variable(data.model, Irji[i=1:data.nbr], start=data.Irji[i]);                     
            data._Iiji=JuMP.@variable(data.model, Iiji[i=1:data.nbr], start=data.Iiji[i]);
        end

        #Adding variables: modeling linearization using piecewise linear interpolation for different ACOPF formulations
        if obj_type == :linear
            JuMP.@variable(data.model, cg[i=1:data.ngen], start=data.cg[i]);
        end

        #Adding variables: modeling shunt FACTS devices as continious switchable shunt succeptance
        if FACTS_bShunt
            data._bss=JuMP.@variable(data.model, bss[i=1:data.nbss], start=data.Bss[i]);
        end

    else
        #Adding variables: modeling voltages for different ACOPF formulations 
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel] 
            data._θ=JuMP.@variable(data.model, θ[1:data.nbus]);                       
            data._V=JuMP.@variable(data.model, V[1:data.nbus]);
        elseif formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel, PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            data._vi=JuMP.@variable(data.model, vi[1:data.nbus]);                       
            data._vr=JuMP.@variable(data.model, vr[1:data.nbus]);
        end
        if formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            data._Wd=JuMP.@variable(data.model, Wd[i=1:data.nbus]);        
        end
        if formulation ∈ [PBRAWVmodel, PNRAWVmodel] 
            data._Wr=JuMP.@variable(data.model, Wr[i=1:data.nbr]);        
            data._Wi=JuMP.@variable(data.model, Wi[i=1:data.nbr]);
        end

        #Adding variables: modeling generator injections for different ACOPF formulations
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
                            PBRARVmodel, CBRARVmodel, PNRARVmodel, 
                            PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            data._pg=JuMP.@variable(data.model, pg[i=1:data.ngen]);        
            data._qg=JuMP.@variable(data.model, qg[i=1:data.ngen]);
        end

        #Adding variables: modeling branch flows for different ACOPF formulations
        if formulation ∈ [PBRAPVmodel, PBRARVmodel, PBRAWVmodel] 
            data._Pij=JuMP.@variable(data.model, Pij[1:data.nbr]);                      
            data._Qij=JuMP.@variable(data.model, Qij[1:data.nbr]);
            data._Pji=JuMP.@variable(data.model, Pji[1:data.nbr]);                     
            data._Qji=JuMP.@variable(data.model, Qji[1:data.nbr]);
        end
        if formulation ∈ [CBRARVmodel, CBRAWVmodel] 
            data._Irij=JuMP.@variable(data.model, Irij[1:data.nbr]);                      
            data._Iiij=JuMP.@variable(data.model, Iiij[1:data.nbr]);
            data._Irji=JuMP.@variable(data.model, Irji[1:data.nbr]);                     
            data._Iiji=JuMP.@variable(data.model, Iiji[1:data.nbr]);
        end

        #Adding variables: modeling linearization using piecewise linear interpolation for different ACOPF formulations
        if obj_type == :linear
            JuMP.@variable(data.model, cg[i=1:data.ngen]);
        end

        #Adding variables: modeling shunt FACTS devices as continious switchable shunt succeptance
        if FACTS_bShunt
            data._bss=JuMP.@variable(data.model, bss[i=1:data.nbss]);
        end
    end

    #Defining objective function
    if obj_type == :linear
        JuMP.@variable(data.model, tgh[i = 1:data.ngen, j = 1:data.lps[i]] >= 0);
        JuMP.@objective(data.model, Min, sum(cg));
        for i=1:data.ngen
            JuMP.@constraint(data.model, cg[i] == sum(data.cgh[i][a] * tgh[i,a] for a = 1:data.lps[i]));
            JuMP.@constraint(data.model, pg[i] == sum(data.pgh[i][a] * tgh[i,a] for a = 1:data.lps[i]));
            JuMP.@constraint(data.model, sum(tgh[i,a] for a=1:data.lps[i]) == 1);
        end
    elseif obj_type == :quadratic
        c0 = sum(data.c0);
    	if data.cost_order == 3
    		JuMP.@objective(data.model, Min, data.c2' * (pg .* pg) + (data.c1)' * pg + c0);
    	elseif data.cost_order == 2
    		JuMP.@objective(data.model, Min, (data.c1)' * pg + c0);
        else
            JuMP.@objective(data.model, Min, c0);
        end
    elseif obj_type == :feasible
        JuMP.@objective(data.model, Min, 0.0);
    end

    #Adding constraints: defining generator injections bounds for different ACOPF formulations
    if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
                        PBRARVmodel, CBRARVmodel, PNRARVmodel, 
                        PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
        JuMP.@constraint(data.model, data.Pmin .<= pg .<= data.Pmax);   
        JuMP.@constraint(data.model, data.Qmin .<= qg .<= data.Qmax);
    end

    #Adding constraint: defining bounds for shunt FACTS devices modeled as continious switchable shunt succeptance
    if FACTS_bShunt
        JuMP.@constraint(data.model, data.bmin .<= bss .<= data.bmax);   
    end

    #Adding constraints: defining voltages bounds for different ACOPF formulations 
    if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel] 
        JuMP.@constraint(data.model, data.Vmin .<= V .<= data.Vmax);   
    elseif formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel, PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
        JuMP.@constraint(data.model, data.Vmin .* data.Vmin .<= vr .* vr + vi .* vi .<= data.Vmax .* data.Vmax);
    end

    if formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
        JuMP.@constraint(data.model, data.Vmin .* data.Vmin .<= Wd .<= data.Vmax .* data.Vmax);
        JuMP.@constraint(data.model, Wd .== vr .* vr + vi .* vi); 
    end
    if formulation ∈ [PBRAWVmodel, PNRAWVmodel] 
        JuMP.@constraint(data.model, Wd .== vr .* vr + vi .* vi);       
        JuMP.@constraint(data.model, Wd .== vr .* vr + vi .* vi);
    end

    for i=1:data.nbus
        PNLexpr = (formulation ∉ [PBRAWVmodel, PNRAWVmodel]) ? JuMP.@NLexpression(data.model, 0.0) : 0.0 ;            
        QNLexpr = (formulation ∉ [PBRAWVmodel, PNRAWVmodel]) ? JuMP.@NLexpression(data.model, 0.0) : 0.0 ;  
        nbal = data.node[string(i)];
            
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
            PBRARVmodel, CBRARVmodel, PNRARVmodel, CBRAWVmodel]
            for g in nbal["gen"]
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr + pg[g]);
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + qg[g]);
            end
        elseif formulation ∈ [PBRAWVmodel, PNRAWVmodel] 
            for g in nbal["gen"]
                PNLexpr += pg[g];
                QNLexpr += qg[g];
            end
        end

        if formulation ∈ [PBRAPVmodel, PBRARVmodel] 
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];         
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr + Pij[k]);
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + Qij[k]);
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr + Pji[k]);
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + Qji[k]);
            end
        elseif formulation ∈ [PBRAWVmodel] 
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];         
                PNLexpr += Pij[k];
                QNLexpr += Qij[k];
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr += Pji[k];
                QNLexpr += Qji[k];
            end
        elseif formulation ∈ [CBRARVmodel, CBRAWVmodel] 
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr + vr[f] * Irij[k] + vi[f] * Iiij[k]);
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + vi[f] * Irij[k] - vr[f] * Iiij[k]);
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr + vr[t] * Irji[k] + vi[t] * Iiji[k]);
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + vi[t] * Irji[k] - vr[t] * Iiji[k]);
            end
        elseif formulation == PNPAPVmodel
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];         
                Y = abs(data.Gf[k] + (data.Bf[k])im);         θY = angle(data.Gf[k] + (data.Bf[k])im); 
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr - data.gf[k] * V[f]^2 - Y * cos(θ[f] - θ[t] - θY) * V[f] * V[t]);
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + data.bf[k] * V[f]^2 - Y * sin(θ[f] - θ[t] - θY) * V[f] * V[t]);
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                Y = abs(data.Gt[k] + (data.Bt[k])im);         θY = angle(data.Gt[k] + (data.Bt[k])im); 
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr - data.gt[k] * V[t]^2 - Y * cos(θ[t] - θ[f] - θY) * V[f] * V[t]);
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + data.bt[k] * V[t]^2 - Y * sin(θ[t] - θ[f] - θY) * V[f] * V[t]);
            end
        elseif formulation == PNRAPVmodel
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr - data.gf[k] * V[f]^2 - (data.Gf[k] * cos(θ[f] - θ[t]) + data.Bf[k] * sin(θ[f] - θ[t])) * V[f] * V[t]);
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + data.bf[k] * V[f]^2 + (data.Bf[k] * cos(θ[f] - θ[t]) - data.Gf[k] * sin(θ[f] - θ[t])) * V[f] * V[t]);
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr - data.gt[k] * V[t]^2 - (data.Gt[k] * cos(θ[t] - θ[f]) + data.Bt[k] * sin(θ[t] - θ[f])) * V[f] * V[t]);
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + data.bt[k] * V[t]^2 + (data.Bt[k] * cos(θ[t] - θ[f]) - data.Gt[k] * sin(θ[t] - θ[f])) * V[f] * V[t]);
            end
        elseif formulation == PNRARVmodel
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];         
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr - data.gf[k] * (vr[f]^2+vi[f]^2) - data.Gf[k] * (vr[f]*vr[t]+vi[f]*vi[t]) + data.Bf[k] * (vr[f]*vi[t]-vi[f]*vr[t]));
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + data.bf[k] * (vr[f]^2+vi[f]^2) + data.Bf[k] * (vr[f]*vr[t]+vi[f]*vi[t]) + data.Gf[k] * (vr[f]*vi[t]-vi[f]*vr[t]));
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr - data.gt[k] * (vr[t]^2+vi[t]^2) - data.Gt[k] * (vr[t]*vr[f]+vi[t]*vi[f]) + data.Bt[k] * (vr[t]*vi[f]-vi[t]*vr[f]));
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + data.bt[k] * (vr[t]^2+vi[t]^2) + data.Bt[k] * (vr[t]*vr[f]+vi[t]*vi[f]) + data.Gt[k] * (vr[t]*vi[f]-vi[t]*vr[f]));
            end
        elseif formulation == PNRAWVmodel
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];         
                PNLexpr += - data.gf[k] * Wd[f] - data.Gf[k] * Wr[k] + data.Bf[k] * Wi[k];
                QNLexpr += + data.bf[k] * Wd[f] + data.Bf[k] * Wr[k] + data.Gf[k] * Wi[k];
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr += - data.gt[k] * Wd[t] - data.Gt[k] * Wr[k] - data.Bt[k] * Wi[k];
                QNLexpr += + data.bt[k] * Wd[t] + data.Bt[k] * Wr[k] - data.Gt[k] * Wi[k];
            end
        end


        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel] 
            PNLexpr = JuMP.@NLexpression(data.model, PNLexpr - data.Gl[i] * V[i]^2);
            QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + data.Bl[i] * V[i]^2);
            if FACTS_bShunt
                for sh in nbal["shunt"]
                    QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + bss[sh] * V[i]^2);
                end
            end
        elseif formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel] 
            PNLexpr = JuMP.@NLexpression(data.model, PNLexpr - data.Gl[i] * (vr[i]^2 + vi[i]^2));
            QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + data.Bl[i] * (vr[i]^2 + vi[i]^2));
            if FACTS_bShunt
                for sh in nbal["shunt"]
                    QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + bss[sh] * (vr[i]^2 + vi[i]^2));
                end
            end
        elseif formulation ∈ [CBRAWVmodel] 
            PNLexpr = JuMP.@NLexpression(data.model, PNLexpr - data.Gl[i] * Wd[i]);
            QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + data.Bl[i] * Wd[i]);
            if FACTS_bShunt
                for sh in nbal["shunt"]
                    QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + bss[sh] * Wd[i]);
                end
            end
        elseif formulation ∈ [PBRAWVmodel, PNRAWVmodel] 
            PNLexpr -= data.Gl[i] * Wd[i];
            QNLexpr += data.Bl[i] * Wd[i];
            if FACTS_bShunt
                for sh in nbal["shunt"]
                    QNLexpr += bss[sh] * Wd[i];
                end
            end
        end

        if formulation ∉ [PBRAWVmodel, PNRAWVmodel]
            JuMP.@NLconstraint(data.model, PNLexpr == data.Pd[i]);
            JuMP.@NLconstraint(data.model, QNLexpr == data.Qd[i]);
        else
            JuMP.@constraint(data.model, PNLexpr == data.Pd[i]);
            JuMP.@constraint(data.model, QNLexpr == data.Qd[i]);
        end
    end
    
    if data.source_type == "matpower" 
        for k=1:data.nbr                   f = data.br[k][1];                         t = data.br[k][2];
            if formulation == PBRAPVmodel
                JuMP.@NLconstraint(data.model, Pij[k] == - data.gf[k] * V[f]^2 - (data.Gf[k] * cos(θ[f] - θ[t]) + data.Bf[k] * sin(θ[f] - θ[t])) * V[f] * V[t]);
                JuMP.@NLconstraint(data.model, Pji[k] == - data.gt[k] * V[t]^2 - (data.Gt[k] * cos(θ[t] - θ[f]) + data.Bt[k] * sin(θ[t] - θ[f])) * V[f] * V[t]);
                JuMP.@NLconstraint(data.model, Qij[k] ==   data.bf[k] * V[f]^2 + (data.Bf[k] * cos(θ[f] - θ[t]) - data.Gf[k] * sin(θ[f] - θ[t])) * V[f] * V[t]); 
                JuMP.@NLconstraint(data.model, Qji[k] ==   data.bt[k] * V[t]^2 + (data.Bt[k] * cos(θ[t] - θ[f]) - data.Gt[k] * sin(θ[t] - θ[f])) * V[f] * V[t]); 
                JuMP.@NLconstraint(data.model, Pij[k]^2 + Qij[k]^2 <= data.Imax[k]^2);
                JuMP.@NLconstraint(data.model, Pji[k]^2 + Qji[k]^2 <= data.Imax[k]^2);
            elseif formulation == PBRARVmodel
                JuMP.@NLconstraint(data.model, Pij[k] == - data.gf[k] * (vr[f]^2+vi[f]^2) - data.Gf[k] * (vr[f]*vr[t]+vi[f]*vi[t]) + data.Bf[k] * (vr[f]*vi[t]-vi[f]*vr[t]));
                JuMP.@NLconstraint(data.model, Pji[k] == - data.gt[k] * (vr[t]^2+vi[t]^2) - data.Gt[k] * (vr[t]*vr[f]+vi[t]*vi[f]) + data.Bt[k] * (vr[t]*vi[f]-vi[t]*vr[f]));
                JuMP.@NLconstraint(data.model, Qij[k] == + data.bf[k] * (vr[f]^2+vi[f]^2) + data.Bf[k] * (vr[f]*vr[t]+vi[f]*vi[t]) + data.Gf[k] * (vr[f]*vi[t]-vi[f]*vr[t])); 
                JuMP.@NLconstraint(data.model, Qji[k] == + data.bt[k] * (vr[t]^2+vi[t]^2) + data.Bt[k] * (vr[t]*vr[f]+vi[t]*vi[f]) + data.Gt[k] * (vr[t]*vi[f]-vi[t]*vr[f])); 
                JuMP.@NLconstraint(data.model, Pij[k]^2 + Qij[k]^2 <= data.Imax[k]^2);
                JuMP.@NLconstraint(data.model, Pji[k]^2 + Qji[k]^2 <= data.Imax[k]^2);
            elseif formulation == PBRAWVmodel
                JuMP.@NLconstraint(data.model, Wr[k] * Wr[k]  + Wi[k] * Wi[k] == Wd[f] * Wd[t]);
                JuMP.@constraint(data.model, Pij[k] == - data.gf[k] * Wd[f] - data.Gf[k] * Wr[k] + data.Bf[k] * Wi[k]);
                JuMP.@constraint(data.model, Pji[k] == - data.gt[k] * Wd[t] - data.Gt[k] * Wr[k] - data.Bt[k] * Wi[k]);
                JuMP.@constraint(data.model, Qij[k] == + data.bf[k] * Wd[f] + data.Bf[k] * Wr[k] + data.Gf[k] * Wi[k]); 
                JuMP.@constraint(data.model, Qji[k] == + data.bt[k] * Wd[t] + data.Bt[k] * Wr[k] - data.Gt[k] * Wi[k]); 
                JuMP.@NLconstraint(data.model, Pij[k]^2 + Qij[k]^2 <= data.Imax[k]^2);
                JuMP.@NLconstraint(data.model, Pji[k]^2 + Qji[k]^2 <= data.Imax[k]^2);
            elseif formulation == CBRARVmodel
                JuMP.@constraint(data.model, Irij[k] == - data.gf[k] * vr[f] + data.bf[k] * vi[f] - data.Gf[k] * vr[t] + data.Bf[k] * vi[t]);
                JuMP.@constraint(data.model, Irji[k] == - data.gt[k] * vr[t] + data.bt[k] * vi[t] - data.Gt[k] * vr[f] + data.Bt[k] * vi[f]);
                JuMP.@constraint(data.model, Iiij[k] == - data.gf[k] * vi[f] - data.bf[k] * vr[f] - data.Gf[k] * vi[t] - data.Bf[k] * vr[t]);
                JuMP.@constraint(data.model, Iiji[k] == - data.gt[k] * vi[t] - data.bt[k] * vr[t] - data.Gt[k] * vi[f] - data.Bt[k] * vr[f]);
                JuMP.@NLconstraint(data.model, Irij[k]^2 * (vr[f]^2+vi[f]^2) + Iiij[k]^2 * (vr[f]^2+vi[f]^2) <= data.Imax[k]^2);
                JuMP.@NLconstraint(data.model, Irji[k]^2 * (vr[t]^2+vi[t]^2) + Iiji[k]^2 * (vr[t]^2+vi[t]^2) <= data.Imax[k]^2); 
            elseif formulation == CBRAWVmodel
                JuMP.@constraint(data.model, Irij[k] == - data.gf[k] * vr[f] + data.bf[k] * vi[f] - data.Gf[k] * vr[t] + data.Bf[k] * vi[t]);
                JuMP.@constraint(data.model, Irji[k] == - data.gt[k] * vr[t] + data.bt[k] * vi[t] - data.Gt[k] * vr[f] + data.Bt[k] * vi[f]);
                JuMP.@constraint(data.model, Iiij[k] == - data.gf[k] * vi[f] - data.bf[k] * vr[f] - data.Gf[k] * vi[t] - data.Bf[k] * vr[t]);
                JuMP.@constraint(data.model, Iiji[k] == - data.gt[k] * vi[t] - data.bt[k] * vr[t] - data.Gt[k] * vi[f] - data.Bt[k] * vr[f]);
                JuMP.@NLconstraint(data.model, Irij[k]^2 * Wd[f] + Iiij[k]^2 * Wd[f] <= data.Imax[k]^2);
                JuMP.@NLconstraint(data.model, Irji[k]^2 * Wd[t] + Iiji[k]^2 * Wd[t] <= data.Imax[k]^2); 
            elseif formulation ∈ [PNPAPVmodel, PNRAPVmodel]
                y¹ = abs2(data.gf[k] + (data.bf[k])im);       y² = abs2(data.Gf[k] + (data.Bf[k])im);       y³ = 2 * abs(data.gf[k] + (data.bf[k])im) * abs(data.Gf[k] + (data.Bf[k])im);
                θ¹ = angle(data.gf[k] + (data.bf[k])im);      θ² = angle(data.Gf[k] + (data.Bf[k])im); 
                JuMP.@NLconstraint(data.model, y¹ * V[f]^4 + y² * V[t]^2 * V[f]^2 + y³ * V[f]^3 * V[t] * cos(θ[f] - θ[t] + θ¹ - θ²) <= data.Imax[k]^2);
                y¹ = abs2(data.gt[k] + (data.bt[k])im);       y² = abs2(data.Gt[k] + (data.Bt[k])im);       y³ = 2 * abs(data.gt[k] + (data.bt[k])im) * abs(data.Gt[k] + (data.Bt[k])im);
                θ¹ = angle(data.gt[k] + (data.bt[k])im);      θ² = angle(data.Gt[k] + (data.Bt[k])im); 
                JuMP.@NLconstraint(data.model, y¹ * V[t]^4 + y² * V[f]^2 * V[t]^2 + y³ * V[t]^3 * V[f] * cos(θ[t] - θ[f] + θ¹ - θ²) <= data.Imax[k]^2);
            elseif formulation == PNRARVmodel
                y = abs2(data.gf[k] + (data.bf[k])im);                    Y = abs2(data.Gf[k] + (data.Bf[k])im);       
                y¹ = 2 * (data.gf[k] * data.Gf[k] + data.bf[k] * data.Bf[k]);   y² = 2 * (data.bf[k] * data.Gf[k] - data.gf[k] * data.Bf[k]); 
                JuMP.@NLconstraint(data.model, y¹ * (vr[f]*vr[t]+vi[f]*vi[t]) * (vr[f]^2+vi[f]^2) + y² * (vr[f]*vi[t]-vi[f]*vr[t]) * (vr[f]^2+vi[f]^2) 
                                    + y * (vr[f]^2+vi[f]^2) * (vr[f]^2+vi[f]^2) + Y * (vr[t]^2+vi[t]^2) * (vr[f]^2+vi[f]^2) <= data.Imax[k]^2);
                y = abs2(data.gt[k] + (data.bt[k])im);                    Y = abs2(data.Gt[k] + (data.Bt[k])im);       
                y¹ = 2 * (data.gt[k] * data.Gt[k] + data.bt[k] * data.Bt[k]);   y² = 2 * (data.bt[k] * data.Gt[k] - data.gt[k] * data.Bt[k]); 
                JuMP.@NLconstraint(data.model, y¹ * (vr[t]*vr[f]+vi[t]*vi[f]) * (vr[t]^2+vi[t]^2) + y² * (vr[t]*vi[f]-vi[t]*vr[f]) * (vr[t]^2+vi[t]^2) 
                                    + y * (vr[t]^2+vi[t]^2) * (vr[t]^2+vi[t]^2) + Y * (vr[f]^2+vi[f]^2) * (vr[t]^2+vi[t]^2) <= data.Imax[k]^2);   
            elseif formulation == PNRAWVmodel
                JuMP.@NLconstraint(data.model, Wr[k] * Wr[k]  + Wi[k] * Wi[k] == Wd[f] * Wd[t]);
                y = abs2(data.gf[k] + (data.bf[k])im);                    Y = abs2(data.Gf[k] + (data.Bf[k])im);       
                y¹ = 2 * (data.gf[k] * data.Gf[k] + data.bf[k] * data.Bf[k]);   y² = 2 * (data.bf[k] * data.Gf[k] - data.gf[k] * data.Bf[k]); 
                JuMP.@constraint(data.model, y¹ * Wr[k] * Wd[f] + y² * Wi[k] * Wd[f] + y * Wd[f] * Wd[f] + Y * Wd[t] * Wd[f] <= data.Imax[k]^2);
                y = abs2(data.gt[k] + (data.bt[k])im);                    Y = abs2(data.Gt[k] + (data.Bt[k])im);       
                y¹ = 2 * (data.gt[k] * data.Gt[k] + data.bt[k] * data.Bt[k]);   y² = 2 * (data.bt[k] * data.Gt[k] - data.gt[k] * data.Bt[k]); 
                JuMP.@constraint(data.model, y¹ * Wr[k] * Wd[t] - y² * Wi[k] * Wd[t] + y * Wd[t] * Wd[t] + Y * Wd[f] * Wd[t] <= data.Imax[k]^2);   
            end
        end
    else    
        for k=1:data.nbr                   f = data.br[k][1];                         t = data.br[k][2];
            if formulation == PBRAPVmodel
                JuMP.@NLconstraint(data.model, Pij[k] == - data.gf[k] * V[f]^2 - (data.Gf[k] * cos(θ[f] - θ[t]) + data.Bf[k] * sin(θ[f] - θ[t])) * V[f] * V[t]);
                JuMP.@NLconstraint(data.model, Pji[k] == - data.gt[k] * V[t]^2 - (data.Gt[k] * cos(θ[t] - θ[f]) + data.Bt[k] * sin(θ[t] - θ[f])) * V[f] * V[t]);
                JuMP.@NLconstraint(data.model, Qij[k] ==   data.bf[k] * V[f]^2 + (data.Bf[k] * cos(θ[f] - θ[t]) - data.Gf[k] * sin(θ[f] - θ[t])) * V[f] * V[t]); 
                JuMP.@NLconstraint(data.model, Qji[k] ==   data.bt[k] * V[t]^2 + (data.Bt[k] * cos(θ[t] - θ[f]) - data.Gt[k] * sin(θ[t] - θ[f])) * V[f] * V[t]); 
                JuMP.@NLconstraint(data.model, Pij[k]^2 + Qij[k]^2 <= data.Imax[k]^2 * V[f]^2);
                JuMP.@NLconstraint(data.model, Pji[k]^2 + Qji[k]^2 <= data.Imax[k]^2 * V[t]^2);
            elseif formulation == PBRARVmodel
                JuMP.@NLconstraint(data.model, Pij[k] == - data.gf[k] * (vr[f]^2+vi[f]^2) - data.Gf[k] * (vr[f]*vr[t]+vi[f]*vi[t]) + data.Bf[k] * (vr[f]*vi[t]-vi[f]*vr[t]));
                JuMP.@NLconstraint(data.model, Pji[k] == - data.gt[k] * (vr[t]^2+vi[t]^2) - data.Gt[k] * (vr[t]*vr[f]+vi[t]*vi[f]) + data.Bt[k] * (vr[t]*vi[f]-vi[t]*vr[f]));
                JuMP.@NLconstraint(data.model, Qij[k] == + data.bf[k] * (vr[f]^2+vi[f]^2) + data.Bf[k] * (vr[f]*vr[t]+vi[f]*vi[t]) + data.Gf[k] * (vr[f]*vi[t]-vi[f]*vr[t])); 
                JuMP.@NLconstraint(data.model, Qji[k] == + data.bt[k] * (vr[t]^2+vi[t]^2) + data.Bt[k] * (vr[t]*vr[f]+vi[t]*vi[f]) + data.Gt[k] * (vr[t]*vi[f]-vi[t]*vr[f])); 
                JuMP.@NLconstraint(data.model, Pij[k]^2 + Qij[k]^2 <= data.Imax[k]^2 * (vr[f]^2+vi[f]^2));
                JuMP.@NLconstraint(data.model, Pji[k]^2 + Qji[k]^2 <= data.Imax[k]^2 * (vr[t]^2+vi[t]^2));
            elseif formulation == PBRAWVmodel
                JuMP.@NLconstraint(data.model, Wr[k] * Wr[k]  + Wi[k] * Wi[k] == Wd[f] * Wd[t]);
                JuMP.@constraint(data.model, Pij[k] == - data.gf[k] * Wd[f] - data.Gf[k] * Wr[k] + data.Bf[k] * Wi[k]);
                JuMP.@constraint(data.model, Pji[k] == - data.gt[k] * Wd[t] - data.Gt[k] * Wr[k] - data.Bt[k] * Wi[k]);
                JuMP.@constraint(data.model, Qij[k] == + data.bf[k] * Wd[f] + data.Bf[k] * Wr[k] + data.Gf[k] * Wi[k]); 
                JuMP.@constraint(data.model, Qji[k] == + data.bt[k] * Wd[t] + data.Bt[k] * Wr[k] - data.Gt[k] * Wi[k]); 
                JuMP.@NLconstraint(data.model, Pij[k]^2 + Qij[k]^2 <= data.Imax[k]^2 * Wd[f]);
                JuMP.@NLconstraint(data.model, Pji[k]^2 + Qji[k]^2 <= data.Imax[k]^2 * Wd[t]);
            elseif formulation == CBRARVmodel
                JuMP.@constraint(data.model, Irij[k] == - data.gf[k] * vr[f] + data.bf[k] * vi[f] - data.Gf[k] * vr[t] + data.Bf[k] * vi[t]);
                JuMP.@constraint(data.model, Irji[k] == - data.gt[k] * vr[t] + data.bt[k] * vi[t] - data.Gt[k] * vr[f] + data.Bt[k] * vi[f]);
                JuMP.@constraint(data.model, Iiij[k] == - data.gf[k] * vi[f] - data.bf[k] * vr[f] - data.Gf[k] * vi[t] - data.Bf[k] * vr[t]);
                JuMP.@constraint(data.model, Iiji[k] == - data.gt[k] * vi[t] - data.bt[k] * vr[t] - data.Gt[k] * vi[f] - data.Bt[k] * vr[f]);
                JuMP.@NLconstraint(data.model, Irij[k]^2 + Iiij[k]^2 <= data.Imax[k]^2);
                JuMP.@NLconstraint(data.model, Irji[k]^2 + Iiji[k]^2 <= data.Imax[k]^2); 
            elseif formulation == CBRAWVmodel
                JuMP.@constraint(data.model, Irij[k] == - data.gf[k] * vr[f] + data.bf[k] * vi[f] - data.Gf[k] * vr[t] + data.Bf[k] * vi[t]);
                JuMP.@constraint(data.model, Irji[k] == - data.gt[k] * vr[t] + data.bt[k] * vi[t] - data.Gt[k] * vr[f] + data.Bt[k] * vi[f]);
                JuMP.@constraint(data.model, Iiij[k] == - data.gf[k] * vi[f] - data.bf[k] * vr[f] - data.Gf[k] * vi[t] - data.Bf[k] * vr[t]);
                JuMP.@constraint(data.model, Iiji[k] == - data.gt[k] * vi[t] - data.bt[k] * vr[t] - data.Gt[k] * vi[f] - data.Bt[k] * vr[f]);
                JuMP.@NLconstraint(data.model, Irij[k]^2 + Iiij[k]^2 <= data.Imax[k]^2);
                JuMP.@NLconstraint(data.model, Irji[k]^2 + Iiji[k]^2 <= data.Imax[k]^2); 
            elseif formulation ∈ [PNPAPVmodel, PNRAPVmodel]
                y¹ = abs2(data.gf[k] + (data.bf[k])im);       y² = abs2(data.Gf[k] + (data.Bf[k])im);       y³ = 2 * abs(data.gf[k] + (data.bf[k])im) * abs(data.Gf[k] + (data.Bf[k])im);
                θ¹ = angle(data.gf[k] + (data.bf[k])im);      θ² = angle(data.Gf[k] + (data.Bf[k])im); 
                JuMP.@NLconstraint(data.model, y¹ * V[f]^2 + y² * V[t]^2 + y³ * V[f] * V[t] * cos(θ[f] - θ[t] + θ¹ - θ²) <= data.Imax[k]^2);
                y¹ = abs2(data.gt[k] + (data.bt[k])im);       y² = abs2(data.Gt[k] + (data.Bt[k])im);       y³ = 2 * abs(data.gt[k] + (data.bt[k])im) * abs(data.Gt[k] + (data.Bt[k])im);
                θ¹ = angle(data.gt[k] + (data.bt[k])im);      θ² = angle(data.Gt[k] + (data.Bt[k])im); 
                JuMP.@NLconstraint(data.model, y¹ * V[t]^2 + y² * V[f]^2 + y³ * V[t] * V[f] * cos(θ[t] - θ[f] + θ¹ - θ²) <= data.Imax[k]^2);
            elseif formulation == PNRARVmodel
                y = abs2(data.gf[k] + (data.bf[k])im);                    Y = abs2(data.Gf[k] + (data.Bf[k])im);       
                y¹ = 2 * (data.gf[k] * data.Gf[k] + data.bf[k] * data.Bf[k]);   y² = 2 * (data.bf[k] * data.Gf[k] - data.gf[k] * data.Bf[k]); 
                JuMP.@NLconstraint(data.model, y¹ * (vr[f]*vr[t]+vi[f]*vi[t]) + y² * (vr[f]*vi[t]-vi[f]*vr[t]) + y * (vr[f]^2+vi[f]^2) + Y * (vr[t]^2+vi[t]^2) <= data.Imax[k]^2);
                y = abs2(data.gt[k] + (data.bt[k])im);                    Y = abs2(data.Gt[k] + (data.Bt[k])im);       
                y¹ = 2 * (data.gt[k] * data.Gt[k] + data.bt[k] * data.Bt[k]);   y² = 2 * (data.bt[k] * data.Gt[k] - data.gt[k] * data.Bt[k]); 
                JuMP.@NLconstraint(data.model, y¹ * (vr[t]*vr[f]+vi[t]*vi[f]) + y² * (vr[t]*vi[f]-vi[t]*vr[f]) + y * (vr[t]^2+vi[t]^2) + Y * (vr[f]^2+vi[f]^2) <= data.Imax[k]^2);   
            elseif formulation == PNRAWVmodel
                JuMP.@NLconstraint(data.model, Wr[k] * Wr[k]  + Wi[k] * Wi[k] == Wd[f] * Wd[t]);
                # JuMP.@NLconstraint(data.model, Wr[k] == vr[f] * vr[t] + vi[f] * vi[t]);   
                # JuMP.@NLconstraint(data.model, Wi[k] == vr[f] * vi[t] - vi[f] * vr[t]);
                
                y = abs2(data.gf[k] + (data.bf[k])im);                    Y = abs2(data.Gf[k] + (data.Bf[k])im);       
                y¹ = 2 * (data.gf[k] * data.Gf[k] + data.bf[k] * data.Bf[k]);   y² = 2 * (data.bf[k] * data.Gf[k] - data.gf[k] * data.Bf[k]); 
                JuMP.@constraint(data.model, y¹ * Wr[k] + y² * Wi[k] + y * Wd[f] + Y * Wd[t] <= data.Imax[k]^2);
                y = abs2(data.gt[k] + (data.bt[k])im);                    Y = abs2(data.Gt[k] + (data.Bt[k])im);       
                y¹ = 2 * (data.gt[k] * data.Gt[k] + data.bt[k] * data.Bt[k]);   y² = 2 * (data.bt[k] * data.Gt[k] - data.gt[k] * data.Bt[k]); 
                JuMP.@constraint(data.model, y¹ * Wr[k] - y² * Wi[k] + y * Wd[t] + Y * Wd[f] <= data.Imax[k]^2);   
            end
        end
    end

    if box_constraints
        # add_box_constraints(data, formulation)
        if formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel, PBRAWVmodel, CBRAWVmodel, PNRAWVmodel]
            JuMP.@constraint(data.model, - data.Vmax .<= vi .<= data.Vmax);                       
            JuMP.@constraint(data.model, - data.Vmax .<= vr .<= data.Vmax); 
        end
        if formulation ∈ [PBRAWVmodel, PNRAWVmodel]
            for k=1:data.nbr                   f = data.br[k][1];                         t = data.br[k][2];
                JuMP.@constraint(data.model, - data.Vmax[f] * data.Vmax[t] <= Wr[k] <= data.Vmax[f] * data.Vmax[t]);
                JuMP.@constraint(data.model, - data.Vmax[f] * data.Vmax[t] <= Wi[k] <= data.Vmax[f] * data.Vmax[t]);
            end
        end
        if data.source_type == "matpower"
            if formulation ∈ [PBRAPVmodel, PBRARVmodel, PBRAWVmodel]
                JuMP.@constraint(data.model, - data.Imax .<= Pij .<= data.Imax);
                JuMP.@constraint(data.model, - data.Imax .<= Pji .<= data.Imax);
                JuMP.@constraint(data.model, - data.Imax .<= Qij .<= data.Imax);
                JuMP.@constraint(data.model, - data.Imax .<= Qji .<= data.Imax);
            elseif formulation ∈ [CBRARVmodel, CBRAWVmodel]
                for k=1:data.nbr                   f = data.br[k][1];                         t = data.br[k][2];
                    JuMP.@constraint(data.model, - data.Imax[k] <= data.Vmin[f] * Irij[k] <= data.Imax[k]);
                    JuMP.@constraint(data.model, - data.Imax[k] <= data.Vmin[f] * Iiij[k] <= data.Imax[k]);
                    JuMP.@constraint(data.model, - data.Imax[k] <= data.Vmin[t] * Irji[k] <= data.Imax[k]);
                    JuMP.@constraint(data.model, - data.Imax[k] <= data.Vmin[t] * Iiji[k] <= data.Imax[k]);
                end
            end
        else
            if formulation ∈ [PBRAPVmodel, PBRARVmodel, PBRAWVmodel]
                for k=1:data.nbr                   f = data.br[k][1];                         t = data.br[k][2];
                    JuMP.@constraint(data.model, - data.Imax[k] * data.Vmax[f] <= Pij[k] <= data.Imax[k] * data.Vmax[f]);
                    JuMP.@constraint(data.model, - data.Imax[k] * data.Vmax[t] <= Pji[k] <= data.Imax[k] * data.Vmax[t]);
                    JuMP.@constraint(data.model, - data.Imax[k] * data.Vmax[f] <= Qij[k] <= data.Imax[k] * data.Vmax[f]);
                    JuMP.@constraint(data.model, - data.Imax[k] * data.Vmax[t] <= Qji[k] <= data.Imax[k] * data.Vmax[t]);
                end
            elseif formulation ∈ [CBRARVmodel, CBRAWVmodel]
                JuMP.@constraint(data.model, - data.Imax .<= Irij .<= data.Imax);
                JuMP.@constraint(data.model, - data.Imax .<= Iiij .<= data.Imax);
                JuMP.@constraint(data.model, - data.Imax .<= Irji .<= data.Imax);
                JuMP.@constraint(data.model, - data.Imax .<= Iiji .<= data.Imax);
            end
        end
    end
end

function create_dcopf_model!(data::PowersenseData;
    formulation = PBDCmodel,
    initialize = true,
    FACTS_bSeries = true,
    obj_type::Symbol = :linear)    #obj_type = "linear" or "quadratic"

    "JuMP Model for ACOPF formulations"
    data.model = JuMP.Model();

    #Variable initialization
    if initialize
        #Adding variables: modeling voltages for DCOPF formulations 
        data._θ=JuMP.@variable(data.model, θ[i=1:data.nbus], start=data.θ[i]);                       

        #Adding variables: modeling generator injections for DCOPF formulation
        data._pg=JuMP.@variable(data.model, pg[i=1:data.ngen], start=data.Pg[i]);        

        #Adding variables: modeling branch flows for different ACOPF formulations
        # if formulation ∈ [PBDCmodel, PNDCmodel] 
        #     JuMP.@variable(data.model, Pij[i=1:data.nbr], start=data.Pij[i]);                      
        #     JuMP.@variable(data.model, Pji[i=1:data.nbr], start=data.Pji[i]);                     
        # end
        data._Pij=JuMP.@variable(data.model, Pij[i=1:data.nbr], start=data.Pij[i]);                      
        data._Pji=JuMP.@variable(data.model, Pji[i=1:data.nbr], start=data.Pji[i]);                     

        #Adding variables: modeling linearization using piecewise linear interpolation for different DCOPF formulations
        if obj_type == :linear
            JuMP.@variable(data.model, cg[i=1:data.ngen], start=data.cg[i]);
        end

    else
        #Adding variables: modeling voltages for different ACOPF formulations 
        data._θ=JuMP.@variable(data.model, θ[1:data.nbus]);                       

        #Adding variables: modeling generator injections for different ACOPF formulations
        data._pg=JuMP.@variable(data.model, pg[i=1:data.ngen]);        

        #Adding variables: modeling branch flows for different ACOPF formulations
        # if formulation ∈ [PBDCmodel, PNDCmodel] 
        #     JuMP.@variable(data.model, Pij[1:data.nbr]);                      
        #     JuMP.@variable(data.model, Pji[1:data.nbr]);                     
        # end
        data._Pij=JuMP.@variable(data.model, Pij[1:data.nbr]);                      
        data._Pji=JuMP.@variable(data.model, Pji[1:data.nbr]);   

        #Adding variables: modeling linearization using piecewise linear interpolation for different ACOPF formulations
        if obj_type == :linear
            JuMP.@variable(data.model, cg[i=1:data.ngen]);
        end
    end

    #Defining objective function
    if obj_type == :linear
        JuMP.@variable(data.model, tgh[i = 1:data.ngen, j = 1:data.lps[i]] >= 0);
        JuMP.@objective(data.model, Min, sum(cg));
        for i=1:data.ngen
            JuMP.@constraint(data.model, cg[i] == sum(data.cgh[i][a] * tgh[i,a] for a = 1:data.lps[i]));
            JuMP.@constraint(data.model, pg[i] == sum(data.pgh[i][a] * tgh[i,a] for a = 1:data.lps[i]));
            JuMP.@constraint(data.model, sum(tgh[i,a] for a=1:data.lps[i]) == 1);
        end
    elseif obj_type == :quadratic
        c0 = sum(data.c0);
    	if data.cost_order == 2
    		JuMP.@objective(data.model, Min, data.c2' * (pg .* pg) + (data.c1)' * pg + c0);
    	elseif data.cost_order == 3
    		JuMP.@objective(data.model, Min, (data.c1)' * pg + c0);
        else
            JuMP.@objective(data.model, Min, c0);
        end
    elseif obj_type == :feasible
        JuMP.@objective(data.model, Min, 0.0);
    end

    #Adding constraints: defining generator injections bounds for different ACOPF formulations
    JuMP.@constraint(data.model, data.Pmin .<= pg .<= data.Pmax);   

    JuMP.@constraint(data.model, Pij .== data.gf + data.Gf + data.Bf .* (data.Aij' * θ) - data.Bf .* (data.Aji' * θ));
    JuMP.@constraint(data.model, Pji .== data.gt + data.Gt + data.Bt .* (data.Aji' * θ) - data.Bt .* (data.Aij' * θ));
    JuMP.@constraint(data.model, data.Ag * pg .== data.Aij * Pij + data.Aji * Pji + data.Pd + data.Gl); 

    JuMP.@constraint(data.model, - data.Imax .<= Pij .<= data.Imax);
    JuMP.@constraint(data.model, - data.Imax .<= Pji .<= data.Imax);
end

# function add_box_constraints(data::PowersenseData, formulation::ACOPF_fromulation; voltage = true, flow = true)
#     if formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel, PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] && voltage
#         JuMP.@constraint(data.model, data.Vmin .<= vi .<= data.Vmax);                       
#         JuMP.@constraint(data.model, data.Vmin .<= vr .<= data.Vmax); 
#     end
#     if formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] && voltage
#         for k=1:data.nbr                   f = data.br[k][1];                         t = data.br[k][2];
#             JuMP.@constraint(data.model, - data.Vmax[f] * data.Vmax[t] <= Wr[k] <= data.Vmax[f] * data.Vmax[t]);
#             JuMP.@constraint(data.model, - data.Vmax[f] * data.Vmax[t] <= Wi[k] <= data.Vmax[f] * data.Vmax[t]);
#         end
#     end
#     if formulation ∈ [PBRAPVmodel, PBRARVmodel, PBRAWVmodel] && flow
#         for k=1:data.nbr                   f = data.br[k][1];                         t = data.br[k][2];
#             JuMP.@constraint(data.model, - data.Imax[k] * data.Vmax[f] <= Pij[k] <= data.Imax[k] * data.Vmax[f]);
#             JuMP.@constraint(data.model, - data.Imax[k] * data.Vmax[t] <= Pji[k] <= data.Imax[k] * data.Vmax[t]);
#             JuMP.@constraint(data.model, - data.Imax[k] * data.Vmax[f] <= Qij[k] <= data.Imax[k] * data.Vmax[f]);
#             JuMP.@constraint(data.model, - data.Imax[k] * data.Vmax[t] <= Qji[k] <= data.Imax[k] * data.Vmax[t]);
#         end
#     elseif formulation ∈ [CBRARVmodel, CBRAWVmodel] && flow
#         JuMP.@constraint(data.model, - data.Imax .<= Irij .<= data.Imax);
#         JuMP.@constraint(data.model, - data.Imax .<= Iiij .<= data.Imax);
#         JuMP.@constraint(data.model, - data.Imax .<= Irji .<= data.Imax);
#         JuMP.@constraint(data.model, - data.Imax .<= Iiji .<= data.Imax);
#     end
# end

function run_opf!(data::PowersenseData;
    formulation::Union{ACOPF_fromulation,DCOPF_fromulation} = PBRAPVmodel,
    create_model = true,
    solver = Powersense.Optimizer,
    initialize = true,
    FACTS_bShunt = true,
    FACTS_bSeries = false,
    box_constraints = false,
    solve = true,
    obj_type::Symbol = :linear)    #obj_type = :linear or :quadratic

    if formulation ∈ [PBDCmodel, PNDCmodel, DCmodel]
        run_dcopf!(data, solver=solver, formulation = formulation, initialize = initialize, FACTS_bSeries = FACTS_bSeries, obj_type = obj_type)  
        return
    end

    if create_model
    	create_opf_model!(data, formulation = formulation, initialize = initialize, FACTS_bShunt = FACTS_bShunt, 
    				FACTS_bSeries = FACTS_bSeries, box_constraints = box_constraints, obj_type = obj_type);
    end
    
    if solve
    	JuMP.set_optimizer(data.model, solver)
    	JuMP.optimize!(data.model);
    	if JuMP.has_values(data.model)
    		data.cost = JuMP.objective_value(data.model)
            if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel] 
                data.θ = JuMP.value.(data._θ)                     
                data.V = JuMP.value.(data._V)
            elseif formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel, PBRAWVmodel, CBRAWVmodel, PNRAWVmodel]
                data.vi = JuMP.value.(data._vi)                     
                data.vr = JuMP.value.(data._vr)
            end
            if formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
                data.Wd = JuMP.value.(data._Wd)    
            end
            if formulation ∈ [PBRAWVmodel, PNRAWVmodel] 
                data.Wr = JuMP.value.(data._Wr)       
                data.Wi = JuMP.value.(data._Wi)
            end
    
            #Adding variables: modeling generator injections for different ACOPF formulations
            if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
                                PBRARVmodel, CBRARVmodel, PNRARVmodel, 
                                PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
                data.Pg = JuMP.value.(data._pg)
                data.Qg = JuMP.value.(data._qg)
            end
    
            #Adding variables: modeling branch flows for different ACOPF formulations
            if formulation ∈ [PBRAPVmodel, PBRARVmodel, PBRAWVmodel] 
                data.Pij = JuMP.value.(data._Pij)                   
                data.Qij = JuMP.value.(data._Qij)
                data.Pji = JuMP.value.(data._Pji)                    
                data.Qji = JuMP.value.(data._Qji)
            end
            if formulation ∈ [CBRARVmodel, CBRAWVmodel] 
                data.Irij = JuMP.value.(data._Irij)                     
                data.Iiij = JuMP.value.(data._Iiij)
                data.Irji = JuMP.value.(data._Irji)                    
                data.Iiji = JuMP.value.(data._Iiji)
            end
    
            #Adding variables: modeling shunt FACTS devices as continious switchable shunt succeptance
            if FACTS_bShunt
                data.Bss = JuMP.value.(data._bss)
            end
    	end
    end
end

function run_opf(data::PowersenseData;
    formulation::Union{ACOPF_fromulation,DCOPF_fromulation} = PBRAPVmodel,
    create_model = true,
    solver = Powersense.Optimizer,
    initialize = true,
    FACTS_bShunt = true,
    FACTS_bSeries = false,
    box_constraints = false,
    solve = true,
    obj_type::Symbol = :linear)    #obj_type = :linear or :quadratic

    if formulation ∈ [PBDCmodel, PNDCmodel, DCmodel]
        run_dcopf(data, solver = solver, formulation = formulation, initialize = initialize, FACTS_bSeries = FACTS_bSeries, obj_type = obj_type)  
        return
    end

    if create_model
    	create_opf_model!(data, formulation = formulation, initialize = initialize, FACTS_bShunt = FACTS_bShunt, 
    				FACTS_bSeries = FACTS_bSeries, box_constraints = box_constraints, obj_type = obj_type);
    end
    
    if solve
    	JuMP.set_optimizer(data.model, solver)
    	JuMP.optimize!(data.model);
    	if JuMP.has_values(data.model)
    		data.cost = JuMP.objective_value(data.model)
    	end
    end
end

function run_dcopf!(data::PowersenseData;
    formulation::DCOPF_fromulation = PBDCmodel,
    create_model = true,
    solver = Powersense.Optimizer,
    initialize = true,
    FACTS_bSeries = false,
    solve = true,
    obj_type::Symbol = :linear)    #obj_type = :linear or :quadratic

    if create_model
    	create_dcopf_model!(data, formulation = formulation, initialize = initialize, FACTS_bSeries = FACTS_bSeries, obj_type = obj_type);
    end
    
    if solve
    	JuMP.set_optimizer(data.model, solver)
    	JuMP.optimize!(data.model);
    	if JuMP.has_values(data.model)
    		data.cost = JuMP.objective_value(data.model)
            data.Pg .= JuMP.value.(data._pg)
            data.θ .= JuMP.value.(data._θ)
            data.V .= ones(length(data.V))
            data.Pij .= JuMP.value.(data._Pij)
            data.Pji .= JuMP.value.(data._Pji)
            data.Wd .= data.V .* data.V;
            data.vr .= data.V .* cos.(data.θ) 
            data.vi .= data.V .* sin.(data.θ) 
    	end
    end
end

function run_dcopf(data::PowersenseData;
    formulation::DCOPF_fromulation = PBDCmodel,
    create_model = true,
    solver = Powersense.Optimizer,
    initialize = true,
    FACTS_bSeries = false,
    solve = true,
    obj_type::Symbol = :linear)    #obj_type = :linear or :quadratic

    if create_model
    	create_dcopf_model!(data, formulation = formulation, initialize = initialize, FACTS_bSeries = FACTS_bSeries, obj_type = obj_type);
    end
    
    if solve
    	JuMP.set_optimizer(data.model, solver)
    	JuMP.optimize!(data.model);
    end
end