# Example for using Real-Time OPF algorithm, aka OptOnline 

import PowerModelsONM as ONM
import PowerModelsDistribution as PMD
import InfrastructureModels as IM
import JuMP
import HiGHS, Ipopt

onm_path = joinpath(dirname(pathof(ONM)), "..")

include("src/OnlineOptDynaGrid.jl")

# #  # # # # #  # # #    Parameters     # #  # # # # #  # # # # #  # # # # #  # # # #
# Coefficients for Exploration Signal omega=(2*pi)/(a0*i+a1)
a0 = 1.0
a1 = 3.0

TimePeriods=1; # number of outer time periods
TimeSteps=60; # number of inner time periods
Vmax=1.05; # Upper voltage limit
Vmin=0.95; # Lower volage limit

# Model free control variables
alpha = 400.0; # Control gain
epsilon = 0.001; # Coefficient on Exploration Signal -- ES = epsilon*cos(omega*t)

# Initialize Dict to hold all data needed for optOnline algorithm
optOnlineData = Dict{String, Any}("a0"=>a0, "a1"=>a1, "TimePeriods"=>TimePeriods, "TimeSteps"=>TimeSteps, 
                "Vmax"=>Vmax, "Vmin"=>Vmin, "alpha"=>alpha, "epsilon"=>epsilon)
# #  # # # # #  # # # # #  # # # # #  # # # # #  # # # # #  # # # # #  # # # # #  # # # # #

eng, _ = ONM.parse_network("$onm_path/test/data/ieee13_feeder.dss")
PMD.apply_voltage_bounds!(eng; vm_lb=Vmin, vm_ub=Vmax)

# # Algorithm does not yet work for storage -- disable all storage
for (s,st) in eng["storage"]
    st["status"] = PMD.DISABLED
    st["pex"] = 0.0
    st["discharge_efficiency"] = 100.0
    st["charge_efficiency"] = 100.0
end

# Modify nominal load
for (i, obj) in eng["load"]
    obj["pd_nom"] *= 0.80
    obj["qd_nom"] *= 0.80
end

math = ONM.transform_data_model(eng)

# PowerModels data for "nominal" network
optOnlineData["pm_baseCase"] = ONM.instantiate_onm_model(math, PMD.LPUBFDiagPowerModel, ONM.build_block_mld)

optOnlineData["solver"] = JuMP.optimizer_with_attributes(
	HiGHS.Optimizer,
	"presolve"=>"off",
	"primal_feasibility_tolerance"=>1e-6,
	"dual_feasibility_tolerance"=>1e-6,
	"mip_feasibility_tolerance"=>1e-4,
	"mip_rel_gap"=>1e-4,
	"small_matrix_value"=>1e-12,
	"allow_unbounded_or_infeasible"=>true
)

optOnlineData["NL_solver"] = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "max_iter"=>6000)

optOnlineData["results"] = Dict{String, Dict{String, Any}}()
optOnlineData["case_math"] = Dict{String, Dict{String, Any}}()
optOnlineData["case_math_bounds"] = Dict{String, Dict{String, Any}}()
for t = 1:TimePeriods

    # If modifying network configuration between time steps T, do so here 




    math = ONM.transform_data_model(eng)

     # Modify cost parameters for generation
     for (g,gen) in math["gen"] 
        if haskey(gen, "cost")
            type = split(gen["source_id"], ".")[1]
            if type == "solar"
                gen["ncost"] = 3
                gen["cost"] = [0.0, 0.10, 0.0]
            elseif type == "voltage_source"
                gen["ncost"] = 3
                gen["cost"] = [0.0, 0.50, 0.0]
            else
                gen["ncost"] = 3
                gen["cost"] = [0.0, 0.35, 0.0]
            end 
        end
    end

    # Need a copy of the data that allows voltage violations for PF simulation
    math_copy = deepcopy(math)
    for (b,bus) in math_copy["bus"]
        bus["vmax"] = fill(Vmax+0.10, size(bus["vmax"]))
        bus["vmin"] = fill(Vmin-0.10, size(bus["vmin"]))
    end
    merge!(optOnlineData["case_math"], Dict("$t" => math_copy))
    merge!(optOnlineData["case_math_bounds"], Dict("$t" => math))

    # Get network configuration
    pm = ONM.instantiate_onm_model(optOnlineData["case_math_bounds"]["$t"], PMD.LPUBFDiagPowerModel, ONM.build_block_mld)
    JuMP.set_optimizer(pm.model, optOnlineData["solver"])
    JuMP.optimize!(pm.model)
    result = IM.build_result(pm, JuMP.solve_time(pm.model); solution_processors=ONM._default_solution_processors)
    sol_eng = PMD.transform_solution(result["solution"], math)
    merge!(optOnlineData["results"], Dict("$t" => result))
end 

# This current function creates a Dict of "pd" and "qd" values for each time step -- currently sets all time step values to the base value taken from the data
optOnlineData = OnlineOptDynaGrid.GenerateLoadTimeSeriesMATH(optOnlineData)

# Modify load time series below to cause load changes
for k in 3:TimeSteps
    for (l, load) in optOnlineData["pm_baseCase"].data["load"]
        if l in ["2", "5", "9", "11","19","23"]
            load["timeInterval"]["1"]["timeStep"][string(k)]["pd"] = load["timeInterval"]["1"]["timeStep"][string(k)]["pd"].*0.90
        end 
    end 
end


# This runs the model-free real-time OPF algorithm and returns the measured bus voltages and all generator setpoints, as well as the OPF at the start of each time period, the connected components(CC), and the DERs within each CC
#tracked_data, resultopf, CC, CC_DER = OnlineOptDynaGrid.runOptOnline(pm_baseCase, case_math, case_math_bounds, results, T, N, a0, a1, alpha, epsilon, Vmax, Vmin, NL_solver);
optOnlineResults = OnlineOptDynaGrid.runOptOnline(optOnlineData)