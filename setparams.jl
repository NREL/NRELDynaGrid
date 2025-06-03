#import OnlineOptDynaGrid

include("src/OnlineOptDynaGrid.jl")

a0 = 5.0
a1 = 5.0
T=3; # number of outer time periods
N=60; # number of inner time periods
N_sim=T*N+1;
Vmax=1.05;
Vmin=0.95;

# Model free control variables
alpha = 400.0; # 400 for CASEMATH
epsilon = 0.001; # smaller than 0.0001 is too slow

Ppv_control=Dict()
Qpv_control=Dict()
Ppv_control["up"]=Dict()
Qpv_control["up"]=Dict()
Ppv_control["down"]=Dict()
Qpv_control["down"]=Dict()
Ppv_control["base"]=Dict()
Qpv_control["base"]=Dict()