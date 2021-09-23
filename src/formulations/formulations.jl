include("type.jl")

function run_OPF(data::PowersenseData;
    formulation = PBRAPVmodel,
    solver = Powersense.Optimizer,
    initialize = false,
    FACTS_bShunt = true,
    FACTS_bSeries = true,
    box_constraints = false,
    solve = true,
    obj_type = "linear")    #obj_type = "linear" or "quadratic"

    "JuMP Model for ACOPF formulations"
    model = JuMP.Model(solver);

    #Variable initialization
    if initialize
        #Adding variables: modeling voltages for different ACOPF formulations 
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel] 
            JuMP.@variable(model, θ[1:data.nbus], start=data.θ);                       
            JuMP.@variable(model, V[1:data.nbus], start=data.V);
        elseif formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel, PBRAWVmodel, CBRAWVmodel, PNRAWVmodel]
            JuMP.@variable(model, vi[1:data.nbus], start=data.vr);                       
            JuMP.@variable(model, vr[1:data.nbus], start=data.vr);
        end
        if formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            JuMP.@variable(model, Wd[i=1:data.nbus], start=data.Wd[i]);        
            JuMP.@variable(model, Wr[i=1:data.nbr], start=data.Wr[i]);        
            JuMP.@variable(model, Wi[i=1:data.nbr], start=data.Wi[i]);
        end

        #Adding variables: modeling generator injections for different ACOPF formulations
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
                            PBRARVmodel, CBRARVmodel, PNRARVmodel, 
                            PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            JuMP.@variable(model, pg[i=1:data.ngen], start=data.Pg[i]);        
            JuMP.@variable(model, qg[i=1:data.ngen], start=data.Qg[i]);
        end

        #Adding variables: modeling branch flows for different ACOPF formulations
        if formulation ∈ [PBRAPVmodel, PBRARVmodel, PBRAWVmodel] 
            JuMP.@variable(model, Pij[1:data.nbr], start=data.Pij);                      
            JuMP.@variable(model, Qij[1:data.nbr], start=data.Pij);
            JuMP.@variable(model, Pji[1:data.nbr], start=data.Pij);                     
            JuMP.@variable(model, Qji[1:data.nbr], start=data.Pij);
        end
        if formulation ∈ [CBRARVmodel, CBRAWVmodel] 
            JuMP.@variable(model, Irij[1:data.nbr], start=data.Irij);                      
            JuMP.@variable(model, Iiij[1:data.nbr], start=data.Iiij);
            JuMP.@variable(model, Irji[1:data.nbr], start=data.Irji);                     
            JuMP.@variable(model, Iiji[1:data.nbr], start=data.Iiji);
        end

        #Adding variables: modeling linearization using piecewise linear interpolation for different ACOPF formulations
        if obj_type == "linear"
            JuMP.@variable(model, cg[i=1:data.ngen], start=data.cg[i]);
        end

        #Adding variables: modeling shunt FACTS devices as continious switchable shunt succeptance
        if FACTS_bShunt
            JuMP.@variable(model, bss[i=1:data.nbss], start=data.Bss[i]);
        end

    else
        #Adding variables: modeling voltages for different ACOPF formulations 
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel] 
            JuMP.@variable(model, θ[1:data.nbus]);                       
            JuMP.@variable(model, V[1:data.nbus]);
        elseif formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel, PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            JuMP.@variable(model, vi[1:data.nbus]);                       
            JuMP.@variable(model, vr[1:data.nbus]);
        end
        if formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            JuMP.@variable(model, Wd[i=1:data.nbus]);        
            JuMP.@variable(model, Wr[i=1:data.nbr]);        
            JuMP.@variable(model, Wi[i=1:data.nbr]);
        end

        #Adding variables: modeling generator injections for different ACOPF formulations
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
                            PBRARVmodel, CBRARVmodel, PNRARVmodel, 
                            PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            JuMP.@variable(model, pg[i=1:data.ngen]);        
            JuMP.@variable(model, qg[i=1:data.ngen]);
        end

        #Adding variables: modeling branch flows for different ACOPF formulations
        if formulation ∈ [PBRAPVmodel, PBRARVmodel, PBRAWVmodel] 
            JuMP.@variable(model, Pij[1:data.nbr]);                      
            JuMP.@variable(model, Qij[1:data.nbr]);
            JuMP.@variable(model, Pji[1:data.nbr]);                     
            JuMP.@variable(model, Qji[1:data.nbr]);
        end
        if formulation ∈ [CBRARVmodel, CBRAWVmodel] 
            JuMP.@variable(model, Irij[1:data.nbr]);                      
            JuMP.@variable(model, Iiij[1:data.nbr]);
            JuMP.@variable(model, Irji[1:data.nbr]);                     
            JuMP.@variable(model, Iiji[1:data.nbr]);
        end

        #Adding variables: modeling linearization using piecewise linear interpolation for different ACOPF formulations
        if obj_type == "linear"
            JuMP.@variable(model, cg[i=1:data.ngen]);
        end

        #Adding variables: modeling shunt FACTS devices as continious switchable shunt succeptance
        if FACTS_bShunt
            JuMP.@variable(model, bss[i=1:data.nbss]);
        end
    end

    #Defining objective function
    if obj_type == "linear"
        JuMP.@variable(model, tgh[i = 1:data.ngen, j = 1:data.lps[i]] >= 0);
        JuMP.@objective(model, Min, sum(cg));
        for i=1:data.ngen
		JuMP.@constraint(model, cg[i] == sum(data.cgh[i][a] * tgh[i,a] for a = 1:data.lps[i]));
		JuMP.@constraint(model, pg[i] == sum(data.pgh[i][a] * tgh[i,a] for a = 1:data.lps[i]));
		JuMP.@constraint(model, sum(tgh[i,a] for a=1:data.lps[i]) == 1);
	end
    elseif formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
                            PBRARVmodel, CBRARVmodel, PNRARVmodel, 
                            PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
        c0 = sum(data.c0);
    	if data.cost_order == 2
    		JuMP.@objective(model, Min, data.c2' * (pg .* pg) + (data.c1)' * pg + c0);
    	elseif data.cost_order == 3
    		JuMP.@objective(model, Min, (data.c1)' * pg + c0);
        else
            JuMP.@objective(model, Min, c0);
        end
    end

    #Adding constraints: defining generator injections bounds for different ACOPF formulations
    if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
                        PBRARVmodel, CBRARVmodel, PNRARVmodel, 
                        PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
        JuMP.@constraint(model, data.Pmin .<= pg .<= data.Pmax);   
        JuMP.@constraint(model, data.Qmin .<= qg .<= data.Qmax);
    end

    #Adding constraint: defining bounds for shunt FACTS devices modeled as continious switchable shunt succeptance
    if FACTS_bShunt
        JuMP.@constraint(model, data.bmin .<= bss .<= data.bmax);   
    end

    #Adding constraints: defining voltages bounds for different ACOPF formulations 
    if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel] 
        JuMP.@constraint(model, data.Vmin .<= V .<= data.Vmax);   
    elseif formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel] 
        JuMP.@constraint(model, data.Vmin .* data.Vmin .<= vr .* vr + vi .* vi .<= data.Vmax .* data.Vmax);
    elseif formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
        JuMP.@constraint(model, data.Vmin .* data.Vmin .<= Wd .<= data.Vmax .* data.Vmax);
        JuMP.@constraint(model, Wd .== vr .* vr + vi .* vi); 
    end

    for i=1:data.nbus
        PNLexpr = JuMP.@NLexpression(model, 0.0);            
        QNLexpr = JuMP.@NLexpression(model, 0.0);    
        nbal = data.node[string(i)];

            
        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel, 
            PBRARVmodel, CBRARVmodel, PNRARVmodel, 
            PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            for g in nbal["gen"]
                PNLexpr = JuMP.@NLexpression(model, PNLexpr + pg[g]);
                QNLexpr = JuMP.@NLexpression(model, QNLexpr + qg[g]);
            end
        end


        if formulation ∈ [PBRAPVmodel, PBRARVmodel, PBRAWVmodel] 
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];         
                PNLexpr = JuMP.@NLexpression(model, PNLexpr + Pij[k]);
                QNLexpr = JuMP.@NLexpression(model, QNLexpr + Qij[k]);
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(model, PNLexpr + Pji[k]);
                QNLexpr = JuMP.@NLexpression(model, QNLexpr + Qji[k]);
            end
        elseif formulation ∈ [CBRARVmodel, CBRAWVmodel] 
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(model, PNLexpr + vr[f] * Irij[k] + vi[f] * Iiij[k]);
                QNLexpr = JuMP.@NLexpression(model, QNLexpr + vi[f] * Irij[k] - vr[f] * Iiij[k]);
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(model, PNLexpr + vr[t] * Irji[k] + vi[t] * Iiji[k]);
                QNLexpr = JuMP.@NLexpression(model, QNLexpr + vi[t] * Irji[k] - vr[t] * Iiji[k]);
            end
        elseif formulation == PNPAPVmodel
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];         
                Y = abs(data.Gf[k] + (data.Bf[k])im);         θY = angle(data.Gf[k] + (data.Bf[k])im); 
                PNLexpr = JuMP.@NLexpression(model, PNLexpr - data.gf[k] * V[f]^2 - Y * cos(θ[f] - θ[t] - θY) * V[f] * V[t]);
                QNLexpr = JuMP.@NLexpression(model, QNLexpr + data.bf[k] * V[f]^2 - Y * sin(θ[f] - θ[t] - θY) * V[f] * V[t]);
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                Y = abs(data.Gt[k] + (data.Bt[k])im);         θY = angle(data.Gt[k] + (data.Bt[k])im); 
                PNLexpr = JuMP.@NLexpression(model, PNLexpr - data.gt[k] * V[t]^2 - Y * cos(θ[t] - θ[f] - θY) * V[f] * V[t]);
                QNLexpr = JuMP.@NLexpression(model, QNLexpr + data.bt[k] * V[t]^2 - Y * sin(θ[t] - θ[f] - θY) * V[f] * V[t]);
            end
        elseif formulation == PNRAPVmodel
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(model, PNLexpr - data.gf[k] * V[f]^2 - (data.Gf[k] * cos(θ[f] - θ[t]) + data.Bf[k] * sin(θ[f] - θ[t])) * V[f] * V[t]);
                QNLexpr = JuMP.@NLexpression(model, QNLexpr + data.bf[k] * V[f]^2 + (data.Bf[k] * cos(θ[f] - θ[t]) - data.Gf[k] * sin(θ[f] - θ[t])) * V[f] * V[t]);
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(model, PNLexpr - data.gt[k] * V[t]^2 - (data.Gt[k] * cos(θ[t] - θ[f]) + data.Bt[k] * sin(θ[t] - θ[f])) * V[f] * V[t]);
                QNLexpr = JuMP.@NLexpression(model, QNLexpr + data.bt[k] * V[t]^2 + (data.Bt[k] * cos(θ[t] - θ[f]) - data.Gt[k] * sin(θ[t] - θ[f])) * V[f] * V[t]);
            end
        elseif formulation == PNRARVmodel
            for k in nbal["Sij"]        f = data.br[k][1];                 t = data.br[k][2];         
                PNLexpr = JuMP.@NLexpression(model, PNLexpr - data.gf[k] * (vr[f]^2+vi[f]^2) - data.Gf[k] * (vr[f]*vr[t]+vi[f]*vi[t]) + data.Bf[k] * (vr[f]*vi[t]-vi[f]*vr[t]));
                QNLexpr = JuMP.@NLexpression(model, QNLexpr + data.bf[k] * (vr[f]^2+vi[f]^2) + data.Bf[k] * (vr[f]*vr[t]+vi[f]*vi[t]) + data.Gf[k] * (vr[f]*vi[t]-vi[f]*vr[t]));
            end
            for k in nbal["Sji"]        f = data.br[k][1];                 t = data.br[k][2];
                PNLexpr = JuMP.@NLexpression(model, PNLexpr - data.gt[k] * (vr[t]^2+vi[t]^2) - data.Gt[k] * (vr[t]*vr[f]+vi[t]*vi[f]) + data.Bt[k] * (vr[t]*vi[f]-vi[t]*vr[f]));
                QNLexpr = JuMP.@NLexpression(model, QNLexpr + data.bt[k] * (vr[t]^2+vi[t]^2) + data.Bt[k] * (vr[t]*vr[f]+vi[t]*vi[f]) + data.Gt[k] * (vr[t]*vi[f]-vi[t]*vr[f]));
            end
        end


        if formulation ∈ [PBRAPVmodel, PNPAPVmodel, PNRAPVmodel] 
            PNLexpr = JuMP.@NLexpression(model, PNLexpr - data.Gl[i] * V[i]^2);
            QNLexpr = JuMP.@NLexpression(model, QNLexpr + data.Bl[i] * V[i]^2);
            if FACTS_bShunt
                for sh in nbal["shunt"]
                    QNLexpr = JuMP.@NLexpression(model, QNLexpr + bss[sh] * V[i]^2);
                end
            end
        elseif formulation ∈ [PBRARVmodel, CBRARVmodel, PNRARVmodel] 
            PNLexpr = JuMP.@NLexpression(model, PNLexpr - data.Gl[i] * (vr[i]^2 + vi[i]^2));
            QNLexpr = JuMP.@NLexpression(model, QNLexpr + data.Bl[i] * (vr[i]^2 + vi[i]^2));
            if FACTS_bShunt
                for sh in nbal["shunt"]
                    QNLexpr = JuMP.@NLexpression(model, QNLexpr + bss[sh] * (vr[i]^2 + vi[i]^2));
                end
            end
        elseif formulation ∈ [PBRAWVmodel, CBRAWVmodel, PNRAWVmodel] 
            PNLexpr = JuMP.@NLexpression(model, PNLexpr - data.Gl[i] * Wd[i]);
            QNLexpr = JuMP.@NLexpression(model, QNLexpr + data.Bl[i] * Wd[i]);
            if FACTS_bShunt
                for sh in nbal["shunt"]
                    QNLexpr = JuMP.@NLexpression(model, QNLexpr + bss[sh] * Wd[i]);
                end
            end
        end

        JuMP.@NLconstraint(model, PNLexpr == data.Pd[i]);
        JuMP.@NLconstraint(model, QNLexpr == data.Qd[i]);
    end
    
    for k=1:data.nbr                   f = data.br[k][1];                         t = data.br[k][2];
        if formulation == PBRAPVmodel
            JuMP.@NLconstraint(model, Pij[k] == - data.gf[k] * V[f]^2 - (data.Gf[k] * cos(θ[f] - θ[t]) + data.Bf[k] * sin(θ[f] - θ[t])) * V[f] * V[t]);
            JuMP.@NLconstraint(model, Pji[k] == - data.gt[k] * V[t]^2 - (data.Gt[k] * cos(θ[t] - θ[f]) + data.Bt[k] * sin(θ[t] - θ[f])) * V[f] * V[t]);
            JuMP.@NLconstraint(model, Qij[k] ==   data.bf[k] * V[f]^2 + (data.Bf[k] * cos(θ[f] - θ[t]) - data.Gf[k] * sin(θ[f] - θ[t])) * V[f] * V[t]); 
            JuMP.@NLconstraint(model, Qji[k] ==   data.bt[k] * V[t]^2 + (data.Bt[k] * cos(θ[t] - θ[f]) - data.Gt[k] * sin(θ[t] - θ[f])) * V[f] * V[t]); 
            JuMP.@NLconstraint(model, Pij[k]^2 + Qij[k]^2 <= data.Imax[k]^2 * V[f]^2);
            JuMP.@NLconstraint(model, Pji[k]^2 + Qji[k]^2 <= data.Imax[k]^2 * V[t]^2);
        elseif formulation == PBRARVmodel
            JuMP.@NLconstraint(model, Pij[k] == - data.gf[k] * (vr[f]^2+vi[f]^2) - data.Gf[k] * (vr[f]*vr[t]+vi[f]*vi[t]) + data.Bf[k] * (vr[f]*vi[t]-vi[f]*vr[t]));
            JuMP.@NLconstraint(model, Pji[k] == - data.gt[k] * (vr[t]^2+vi[t]^2) - data.Gt[k] * (vr[t]*vr[f]+vi[t]*vi[f]) + data.Bt[k] * (vr[t]*vi[f]-vi[t]*vr[f]));
            JuMP.@NLconstraint(model, Qij[k] == + data.bf[k] * (vr[f]^2+vi[f]^2) + data.Bf[k] * (vr[f]*vr[t]+vi[f]*vi[t]) + data.Gf[k] * (vr[f]*vi[t]-vi[f]*vr[t])); 
            JuMP.@NLconstraint(model, Qji[k] == + data.bt[k] * (vr[t]^2+vi[t]^2) + data.Bt[k] * (vr[t]*vr[f]+vi[t]*vi[f]) + data.Gt[k] * (vr[t]*vi[f]-vi[t]*vr[f])); 
            JuMP.@NLconstraint(model, Pij[k]^2 + Qij[k]^2 <= data.Imax[k]^2 * (vr[f]^2+vi[f]^2));
            JuMP.@NLconstraint(model, Pji[k]^2 + Qji[k]^2 <= data.Imax[k]^2 * (vr[t]^2+vi[t]^2));
        elseif formulation == CBRARVmodel
            JuMP.@constraint(model, Irij[k] == - data.gf[k] * vr[f] + data.bf[k] * vi[f] - data.Gf[k] * vr[t] + data.Bf[k] * vi[t]);
            JuMP.@constraint(model, Irji[k] == - data.gt[k] * vr[t] + data.bt[k] * vi[t] - data.Gt[k] * vr[f] + data.Bt[k] * vi[f]);
            JuMP.@constraint(model, Iiij[k] == - data.gf[k] * vi[f] - data.bf[k] * vr[f] - data.Gf[k] * vi[t] - data.Bf[k] * vr[t]);
            JuMP.@constraint(model, Iiji[k] == - data.gt[k] * vi[t] - data.bt[k] * vr[t] - data.Gt[k] * vi[f] - data.Bt[k] * vr[f]);
            JuMP.@NLconstraint(model, Irij[k]^2 + Iiij[k]^2 <= data.Imax[k]^2);
            JuMP.@NLconstraint(model, Irji[k]^2 + Iiji[k]^2 <= data.Imax[k]^2); 
        elseif formulation == CBRAWVmodel
            JuMP.@constraint(model, Irij[k] == - data.gf[k] * vr[f] + data.bf[k] * vi[f] - data.Gf[k] * vr[t] + data.Bf[k] * vi[t]);
            JuMP.@constraint(model, Irji[k] == - data.gt[k] * vr[t] + data.bt[k] * vi[t] - data.Gt[k] * vr[f] + data.Bt[k] * vi[f]);
            JuMP.@constraint(model, Iiij[k] == - data.gf[k] * vi[f] - data.bf[k] * vr[f] - data.Gf[k] * vi[t] - data.Bf[k] * vr[t]);
            JuMP.@constraint(model, Iiji[k] == - data.gt[k] * vi[t] - data.bt[k] * vr[t] - data.Gt[k] * vi[f] - data.Bt[k] * vr[f]);
            JuMP.@NLconstraint(model, Irij[k]^2 + Iiij[k]^2 <= data.Imax[k]^2);
            JuMP.@NLconstraint(model, Irji[k]^2 + Iiji[k]^2 <= data.Imax[k]^2); 
        elseif formulation ∈ [PNPAPVmodel, PNRAPVmodel]
            y¹ = abs2(data.gf[k] + (data.bf[k])im);       y² = abs2(data.Gf[k] + (data.Bf[k])im);       y³ = 2 * abs(data.gf[k] + (data.bf[k])im) * abs(data.Gf[k] + (data.Bf[k])im);
            θ¹ = angle(data.gf[k] + (data.bf[k])im);      θ² = angle(data.Gf[k] + (data.Bf[k])im); 
            JuMP.@NLconstraint(model, y¹ * V[f]^2 + y² * V[t]^2 + y³ * V[f] * V[t] * cos(θ[f] - θ[t] + θ¹ - θ²) <= data.Imax[k]^2);
            y¹ = abs2(data.gt[k] + (data.bt[k])im);       y² = abs2(data.Gt[k] + (data.Bt[k])im);       y³ = 2 * abs(data.gt[k] + (data.bt[k])im) * abs(data.Gt[k] + (data.Bt[k])im);
            θ¹ = angle(data.gt[k] + (data.bt[k])im);      θ² = angle(data.Gt[k] + (data.Bt[k])im); 
            JuMP.@NLconstraint(model, y¹ * V[t]^2 + y² * V[f]^2 + y³ * V[t] * V[f] * cos(θ[t] - θ[f] + θ¹ - θ²) <= data.Imax[k]^2);
        elseif formulation == PNRARVmodel
            y = abs2(data.gf[k] + (data.bf[k])im);                    Y = abs2(data.Gf[k] + (data.Bf[k])im);       
            y¹ = 2 * (data.gf[k] * data.Gf[k] + data.bf[k] * data.Bf[k]);   y² = 2 * (data.bf[k] * data.Gf[k] - data.gf[k] * data.Bf[k]); 
            JuMP.@NLconstraint(model, y¹ * (vr[f]*vr[t]+vi[f]*vi[t]) + y² * (vr[f]*vi[t]-vi[f]*vr[t]) + y * (vr[f]^2+vi[f]^2) + Y * (vr[t]^2+vi[t]^2) <= data.Imax[k]^2);
            y = abs2(data.gt[k] + (data.bt[k])im);                    Y = abs2(data.Gt[k] + (data.Bt[k])im);       
            y¹ = 2 * (data.gt[k] * data.Gt[k] + data.bt[k] * data.Bt[k]);   y² = 2 * (data.bt[k] * data.Gt[k] - data.gt[k] * data.Bt[k]); 
            JuMP.@NLconstraint(model, y¹ * (vr[t]*vr[f]+vi[t]*vi[f]) + y² * (vr[t]*vi[f]-vi[t]*vr[f]) + y * (vr[t]^2+vi[t]^2) + Y * (vr[f]^2+vi[f]^2) <= data.Imax[k]^2);   
        end
    end
    
    if solve
    	JuMP.optimize!(model);
    	if has_values(model)
    		data.cost = JuMP.objective_value(model)
    	end
    end
end
