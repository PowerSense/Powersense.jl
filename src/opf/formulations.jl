include("type.jl")

function create_opf_model!(data::PowersenseData;
    formulation = PBRAPVmodel,
    initialize = false,
    FACTS_bShunt = true,
    FACTS_bSeries = true,
    box_constraints = false,
    obj_type = "linear")    #obj_type = "linear" or "quadratic"

    "JuMP Model for ACOPF formulations"
    data.model = JuMP.Model();

    #Variable initialization
    if initialize
        #Adding variables: modeling voltages for different ACOPF formulations 
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel] 
            JuMP.@variable(data.model, θ[1:data.nbus], start=data.θ);                       
            JuMP.@variable(data.model, V[1:data.nbus], start=data.V);
        elseif formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel, PBRAWVmodel, CBRAWVmodel, PNRAWVmodel]
            JuMP.@variable(data.model, vi[1:data.nbus], start=data.vi);                       
            JuMP.@variable(data.model, vr[1:data.nbus], start=data.vr);
        end
        if formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            JuMP.@variable(data.model, Wd[i=1:data.nbus], start=data.Wd[i]);        
            JuMP.@variable(data.model, Wr[i=1:data.nbr], start=data.Wr[i]);        
            JuMP.@variable(data.model, Wi[i=1:data.nbr], start=data.Wi[i]);
        end

        #Adding variables: modeling generator injections for different ACOPF formulations
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
                            PBRARVmodel, CBRARVmodel, PNRARVmodel, 
                            PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            JuMP.@variable(data.model, pg[i=1:data.ngen], start=data.Pg[i]);        
            JuMP.@variable(data.model, qg[i=1:data.ngen], start=data.Qg[i]);
        end

        #Adding variables: modeling branch flows for different ACOPF formulations
        if formulation ∈ [PBRAPVmodel, PBRARVmodel, PBRAWVmodel] 
            JuMP.@variable(data.model, Pij[1:data.nbr], start=data.Pij);                      
            JuMP.@variable(data.model, Qij[1:data.nbr], start=data.Qij);
            JuMP.@variable(data.model, Pji[1:data.nbr], start=data.Pji);                     
            JuMP.@variable(data.model, Qji[1:data.nbr], start=data.Qji);
        end
        if formulation ∈ [CBRARVmodel, CBRAWVmodel] 
            JuMP.@variable(data.model, Irij[1:data.nbr], start=data.Irij);                      
            JuMP.@variable(data.model, Iiij[1:data.nbr], start=data.Iiij);
            JuMP.@variable(data.model, Irji[1:data.nbr], start=data.Irji);                     
            JuMP.@variable(data.model, Iiji[1:data.nbr], start=data.Iiji);
        end

        #Adding variables: modeling linearization using piecewise linear interpolation for different ACOPF formulations
        if obj_type == "linear"
            JuMP.@variable(data.model, cg[i=1:data.ngen], start=data.cg[i]);
        end

        #Adding variables: modeling shunt FACTS devices as continious switchable shunt succeptance
        if FACTS_bShunt
            JuMP.@variable(data.model, bss[i=1:data.nbss], start=data.Bss[i]);
        end

    else
        #Adding variables: modeling voltages for different ACOPF formulations 
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel] 
            JuMP.@variable(data.model, θ[1:data.nbus]);                       
            JuMP.@variable(data.model, V[1:data.nbus]);
        elseif formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel, PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            JuMP.@variable(data.model, vi[1:data.nbus]);                       
            JuMP.@variable(data.model, vr[1:data.nbus]);
        end
        if formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            JuMP.@variable(data.model, Wd[i=1:data.nbus]);        
            JuMP.@variable(data.model, Wr[i=1:data.nbr]);        
            JuMP.@variable(data.model, Wi[i=1:data.nbr]);
        end

        #Adding variables: modeling generator injections for different ACOPF formulations
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
                            PBRARVmodel, CBRARVmodel, PNRARVmodel, 
                            PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            JuMP.@variable(data.model, pg[i=1:data.ngen]);        
            JuMP.@variable(data.model, qg[i=1:data.ngen]);
        end

        #Adding variables: modeling branch flows for different ACOPF formulations
        if formulation ∈ [PBRAPVmodel, PBRARVmodel, PBRAWVmodel] 
            JuMP.@variable(data.model, Pij[1:data.nbr]);                      
            JuMP.@variable(data.model, Qij[1:data.nbr]);
            JuMP.@variable(data.model, Pji[1:data.nbr]);                     
            JuMP.@variable(data.model, Qji[1:data.nbr]);
        end
        if formulation ∈ [CBRARVmodel, CBRAWVmodel] 
            JuMP.@variable(data.model, Irij[1:data.nbr]);                      
            JuMP.@variable(data.model, Iiij[1:data.nbr]);
            JuMP.@variable(data.model, Irji[1:data.nbr]);                     
            JuMP.@variable(data.model, Iiji[1:data.nbr]);
        end

        #Adding variables: modeling linearization using piecewise linear interpolation for different ACOPF formulations
        if obj_type == "linear"
            JuMP.@variable(data.model, cg[i=1:data.ngen]);
        end

        #Adding variables: modeling shunt FACTS devices as continious switchable shunt succeptance
        if FACTS_bShunt
            JuMP.@variable(data.model, bss[i=1:data.nbss]);
        end
    end

    #Defining objective function
    if obj_type == "linear"
        JuMP.@variable(data.model, tgh[i = 1:data.ngen, j = 1:data.lps[i]] >= 0);
        JuMP.@objective(data.model, Min, sum(cg));
        for i=1:data.ngen
            JuMP.@constraint(data.model, cg[i] == sum(data.cgh[i][a] * tgh[i,a] for a = 1:data.lps[i]));
            JuMP.@constraint(data.model, pg[i] == sum(data.pgh[i][a] * tgh[i,a] for a = 1:data.lps[i]));
            JuMP.@constraint(data.model, sum(tgh[i,a] for a=1:data.lps[i]) == 1);
        end
    elseif obj_type == "quadratic"
        c0 = sum(data.c0);
    	if data.cost_order == 2
    		JuMP.@objective(data.model, Min, data.c2' * (pg .* pg) + (data.c1)' * pg + c0);
    	elseif data.cost_order == 3
    		JuMP.@objective(data.model, Min, (data.c1)' * pg + c0);
        else
            JuMP.@objective(data.model, Min, c0);
        end
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

    for i=1:data.nbus
        PNLexpr = JuMP.@NLexpression(data.model, 0.0);            
        QNLexpr = JuMP.@NLexpression(data.model, 0.0);    
        nbal = data.node[string(i)];
            
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
            PBRARVmodel, CBRARVmodel, PNRARVmodel, 
            PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            for g in nbal["gen"]
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr + pg[g]);
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + qg[g]);
            end
        end

        if formulation ∈ [PBRAPVmodel, PBRARVmodel, PBRAWVmodel] 
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];         
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr + Pij[k]);
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + Qij[k]);
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(data.model, PNLexpr + Pji[k]);
                QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + Qji[k]);
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
        elseif formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            PNLexpr = JuMP.@NLexpression(data.model, PNLexpr - data.Gl[i] * Wd[i]);
            QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + data.Bl[i] * Wd[i]);
            if FACTS_bShunt
                for sh in nbal["shunt"]
                    QNLexpr = JuMP.@NLexpression(data.model, QNLexpr + bss[sh] * Wd[i]);
                end
            end
        end

        JuMP.@NLconstraint(data.model, PNLexpr == data.Pd[i]);
        JuMP.@NLconstraint(data.model, QNLexpr == data.Qd[i]);
    end
    
    # for k=1:data.nbr                   f = data.br[k][1];                         t = data.br[k][2];
    #     if formulation == PBRAPVmodel
    #         JuMP.@NLconstraint(data.model, Pij[k] == - data.gf[k] * V[f]^2 - (data.Gf[k] * cos(θ[f] - θ[t]) + data.Bf[k] * sin(θ[f] - θ[t])) * V[f] * V[t]);
    #         JuMP.@NLconstraint(data.model, Pji[k] == - data.gt[k] * V[t]^2 - (data.Gt[k] * cos(θ[t] - θ[f]) + data.Bt[k] * sin(θ[t] - θ[f])) * V[f] * V[t]);
    #         JuMP.@NLconstraint(data.model, Qij[k] ==   data.bf[k] * V[f]^2 + (data.Bf[k] * cos(θ[f] - θ[t]) - data.Gf[k] * sin(θ[f] - θ[t])) * V[f] * V[t]); 
    #         JuMP.@NLconstraint(data.model, Qji[k] ==   data.bt[k] * V[t]^2 + (data.Bt[k] * cos(θ[t] - θ[f]) - data.Gt[k] * sin(θ[t] - θ[f])) * V[f] * V[t]); 
    #         JuMP.@NLconstraint(data.model, Pij[k]^2 + Qij[k]^2 <= data.Imax[k]^2 * V[f]^2);
    #         JuMP.@NLconstraint(data.model, Pji[k]^2 + Qji[k]^2 <= data.Imax[k]^2 * V[t]^2);
    #     elseif formulation == PBRARVmodel
    #         JuMP.@NLconstraint(data.model, Pij[k] == - data.gf[k] * (vr[f]^2+vi[f]^2) - data.Gf[k] * (vr[f]*vr[t]+vi[f]*vi[t]) + data.Bf[k] * (vr[f]*vi[t]-vi[f]*vr[t]));
    #         JuMP.@NLconstraint(data.model, Pji[k] == - data.gt[k] * (vr[t]^2+vi[t]^2) - data.Gt[k] * (vr[t]*vr[f]+vi[t]*vi[f]) + data.Bt[k] * (vr[t]*vi[f]-vi[t]*vr[f]));
    #         JuMP.@NLconstraint(data.model, Qij[k] == + data.bf[k] * (vr[f]^2+vi[f]^2) + data.Bf[k] * (vr[f]*vr[t]+vi[f]*vi[t]) + data.Gf[k] * (vr[f]*vi[t]-vi[f]*vr[t])); 
    #         JuMP.@NLconstraint(data.model, Qji[k] == + data.bt[k] * (vr[t]^2+vi[t]^2) + data.Bt[k] * (vr[t]*vr[f]+vi[t]*vi[f]) + data.Gt[k] * (vr[t]*vi[f]-vi[t]*vr[f])); 
    #         JuMP.@NLconstraint(data.model, Pij[k]^2 + Qij[k]^2 <= data.Imax[k]^2 * (vr[f]^2+vi[f]^2));
    #         JuMP.@NLconstraint(data.model, Pji[k]^2 + Qji[k]^2 <= data.Imax[k]^2 * (vr[t]^2+vi[t]^2));
    #     elseif formulation == CBRARVmodel
    #         JuMP.@constraint(data.model, Irij[k] == - data.gf[k] * vr[f] + data.bf[k] * vi[f] - data.Gf[k] * vr[t] + data.Bf[k] * vi[t]);
    #         JuMP.@constraint(data.model, Irji[k] == - data.gt[k] * vr[t] + data.bt[k] * vi[t] - data.Gt[k] * vr[f] + data.Bt[k] * vi[f]);
    #         JuMP.@constraint(data.model, Iiij[k] == - data.gf[k] * vi[f] - data.bf[k] * vr[f] - data.Gf[k] * vi[t] - data.Bf[k] * vr[t]);
    #         JuMP.@constraint(data.model, Iiji[k] == - data.gt[k] * vi[t] - data.bt[k] * vr[t] - data.Gt[k] * vi[f] - data.Bt[k] * vr[f]);
    #         JuMP.@NLconstraint(data.model, Irij[k]^2 + Iiij[k]^2 <= data.Imax[k]^2);
    #         JuMP.@NLconstraint(data.model, Irji[k]^2 + Iiji[k]^2 <= data.Imax[k]^2); 
    #     elseif formulation == CBRAWVmodel
    #         JuMP.@constraint(data.model, Irij[k] == - data.gf[k] * vr[f] + data.bf[k] * vi[f] - data.Gf[k] * vr[t] + data.Bf[k] * vi[t]);
    #         JuMP.@constraint(data.model, Irji[k] == - data.gt[k] * vr[t] + data.bt[k] * vi[t] - data.Gt[k] * vr[f] + data.Bt[k] * vi[f]);
    #         JuMP.@constraint(data.model, Iiij[k] == - data.gf[k] * vi[f] - data.bf[k] * vr[f] - data.Gf[k] * vi[t] - data.Bf[k] * vr[t]);
    #         JuMP.@constraint(data.model, Iiji[k] == - data.gt[k] * vi[t] - data.bt[k] * vr[t] - data.Gt[k] * vi[f] - data.Bt[k] * vr[f]);
    #         JuMP.@NLconstraint(data.model, Irij[k]^2 + Iiij[k]^2 <= data.Imax[k]^2);
    #         JuMP.@NLconstraint(data.model, Irji[k]^2 + Iiji[k]^2 <= data.Imax[k]^2); 
    #     elseif formulation ∈ [PNPAPVmodel, PNRAPVmodel]
    #         y¹ = abs2(data.gf[k] + (data.bf[k])im);       y² = abs2(data.Gf[k] + (data.Bf[k])im);       y³ = 2 * abs(data.gf[k] + (data.bf[k])im) * abs(data.Gf[k] + (data.Bf[k])im);
    #         θ¹ = angle(data.gf[k] + (data.bf[k])im);      θ² = angle(data.Gf[k] + (data.Bf[k])im); 
    #         JuMP.@NLconstraint(data.model, y¹ * V[f]^2 + y² * V[t]^2 + y³ * V[f] * V[t] * cos(θ[f] - θ[t] + θ¹ - θ²) <= data.Imax[k]^2);
    #         y¹ = abs2(data.gt[k] + (data.bt[k])im);       y² = abs2(data.Gt[k] + (data.Bt[k])im);       y³ = 2 * abs(data.gt[k] + (data.bt[k])im) * abs(data.Gt[k] + (data.Bt[k])im);
    #         θ¹ = angle(data.gt[k] + (data.bt[k])im);      θ² = angle(data.Gt[k] + (data.Bt[k])im); 
    #         JuMP.@NLconstraint(data.model, y¹ * V[t]^2 + y² * V[f]^2 + y³ * V[t] * V[f] * cos(θ[t] - θ[f] + θ¹ - θ²) <= data.Imax[k]^2);
    #     elseif formulation == PNRARVmodel
    #         y = abs2(data.gf[k] + (data.bf[k])im);                    Y = abs2(data.Gf[k] + (data.Bf[k])im);       
    #         y¹ = 2 * (data.gf[k] * data.Gf[k] + data.bf[k] * data.Bf[k]);   y² = 2 * (data.bf[k] * data.Gf[k] - data.gf[k] * data.Bf[k]); 
    #         JuMP.@NLconstraint(data.model, y¹ * (vr[f]*vr[t]+vi[f]*vi[t]) + y² * (vr[f]*vi[t]-vi[f]*vr[t]) + y * (vr[f]^2+vi[f]^2) + Y * (vr[t]^2+vi[t]^2) <= data.Imax[k]^2);
    #         y = abs2(data.gt[k] + (data.bt[k])im);                    Y = abs2(data.Gt[k] + (data.Bt[k])im);       
    #         y¹ = 2 * (data.gt[k] * data.Gt[k] + data.bt[k] * data.Bt[k]);   y² = 2 * (data.bt[k] * data.Gt[k] - data.gt[k] * data.Bt[k]); 
    #         JuMP.@NLconstraint(data.model, y¹ * (vr[t]*vr[f]+vi[t]*vi[f]) + y² * (vr[t]*vi[f]-vi[t]*vr[f]) + y * (vr[t]^2+vi[t]^2) + Y * (vr[f]^2+vi[f]^2) <= data.Imax[k]^2);   
    #     end
    # end

    if box_constraints
        # add_box_constraints(data, formulation)
        if formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel, PBRAWVmodel, CBRAWVmodel, PNRAWVmodel]
            JuMP.@constraint(data.model, data.Vmin .<= vi .<= data.Vmax);                       
            JuMP.@constraint(data.model, data.Vmin .<= vr .<= data.Vmax); 
        end
        if formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel]
            for k=1:data.nbr                   f = data.br[k][1];                         t = data.br[k][2];
                JuMP.@constraint(data.model, - data.Vmax[f] * data.Vmax[t] <= Wr[k] <= data.Vmax[f] * data.Vmax[t]);
                JuMP.@constraint(data.model, - data.Vmax[f] * data.Vmax[t] <= Wi[k] <= data.Vmax[f] * data.Vmax[t]);
            end
        end
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
    formulation = PBRAPVmodel,
    create_model = true,
    solver = Powersense.Optimizer,
    initialize = false,
    FACTS_bShunt = true,
    FACTS_bSeries = false,
    box_constraints = false,
    solve = true,
    obj_type = "linear")    #obj_type = "linear" or "quadratic"

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
