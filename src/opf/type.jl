mutable struct ACOPF_fromulation ACOPF_fromulation() = new() end
mutable struct DCOPF_fromulation DCOPF_fromulation() = new() end

export PBRAPVmodel, PBRARVmodel, PBRAWVmodel, CBRARVmodel, CBRAWVmodel, PNPAPVmodel, PNRAPVmodel, PNRARVmodel, PNRAWVmodel
export DCmodel, PNDCmodel, PBDCmodel

PBRAPVmodel =  ACOPF_fromulation();     #Power Branch Flow Rectangular Admittance Polar Voltage
PBRARVmodel =  ACOPF_fromulation();     #Power Branch Flow Rectangular Admittance Rectangular Voltage
PBRAWVmodel =  ACOPF_fromulation();     #Power Branch Flow Rectangular Admittance W-matrix Voltage
CBRARVmodel =  ACOPF_fromulation();     #Current Branch Flow Rectangular Admittance Rectangular Voltage
CBRAWVmodel =  ACOPF_fromulation();     #Current Branch Flow Rectangular Admittance W-matrix Voltage
PNPAPVmodel =  ACOPF_fromulation();     #Power Nodal Flow Polar Admittance Polar Voltage
PNRAPVmodel =  ACOPF_fromulation();     #Power Nodal Flow Rectangular Admittance Polar Voltage
PNRARVmodel =  ACOPF_fromulation();     #Power Nodal Flow Rectangular Admittance Rectangular Voltage
PNRAWVmodel =  ACOPF_fromulation();     #Power Nodal Flow Rectangular Admittance W-matrix Voltage

PNDCmodel =  DCOPF_fromulation();       #Power Branch-Flow DC Approximation 
PBDCmodel =  DCOPF_fromulation();       #Power Nodal DC Approximation 
DCmodel =  DCOPF_fromulation();         #Power Nodal DC Approximation 
