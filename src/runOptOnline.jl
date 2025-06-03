# Function to run OptOnline Algorithm given a PowerModels model (pm), mathematical data with and without tight voltage bounds (case_math_bounds, case_math), 
#           starting OPF results for each time period (results), the number of time periods (T), and the number of steps within each time period (N)


function runOptOnline(optOnlineData::Dict{String,Any})
    pm_baseCase = optOnlineData["pm_baseCase"];
    T = optOnlineData["TimePeriods"];
    N = optOnlineData["TimeSteps"];
    epsilon = optOnlineData["epsilon"];
    alpha = optOnlineData["alpha"];

    case_math = optOnlineData["case_math"];
    case_math_bounds = optOnlineData["case_math_bounds"];

    optOnlineData["N_sim"]=T*N+1;
    terminals=pm_baseCase.data["conductor_ids"]
    map_id_pairs = Dict{Int,Tuple{Int,Int}}(id => (PMD.ref(pm_baseCase, 0, :bus_block_map, sw["f_bus"]),PMD.ref(pm_baseCase, 0, :bus_block_map, sw["t_bus"])) for (id,sw) in PMD.ref(pm_baseCase, 0, :switch))
    gen_block_map = PMD.ref(pm_baseCase, 0, :gen_block_map)
    sbase = pm_baseCase.data["settings"]["sbase_default"]

    optOnlineData = OnlineOptDynaGrid.fixDataToSolution(optOnlineData)

    optOnlineData = OnlineOptDynaGrid.initializeData(optOnlineData)
    num_vol = optOnlineData["num_vol"]

    optOnlineData["resultopf"] = Dict{String, Any}("$t" => Dict{String, Any}() for t in 1:T)
    for t=1:T
        delete!(case_math_bounds["$t"], "options")
        delete!(case_math["$t"], "options")

        global optOnlineData["resultopf"]["$t"]=PMD.solve_mc_model(case_math_bounds["$t"], PMD.IVRUPowerModel, optOnlineData["NL_solver"], OnlineOptDynaGrid.build_mc_opf_oltc_dg; solution_processors=[PMD.sol_data_model!], make_si=false)
        #global resultopf["$t"]=PMD.solve_mc_model(case_math_bounds["$t"], PMD.ACPUPowerModel, NL_solver, OnlineOptDynaGrid.build_mc_opf_oltc_dg; solution_processors=[PMD.sol_data_model!], make_si=false)
        @show optOnlineData["resultopf"]["$t"]["termination_status"] 
    end

    optOnlineData = OnlineOptDynaGrid.fill_lam_id(optOnlineData)

    optOnlineData = OnlineOptDynaGrid.sampleProbeSignal(optOnlineData)

    for t=1:T

        optOnlineData = OnlineOptDynaGrid.fixGen(optOnlineData, t)

        rt=(t-1)*N+1:t*N   
        for n in rt      
            @show n
            for (l,load) in pm_baseCase.data["load"]
                optOnlineData["tracked_data"]["$t"]["load"][l]["pd"][:,n] = load["timeInterval"]["$t"]["timeStep"][string(n-N*(t-1))]["pd"]*sbase
                optOnlineData["tracked_data"]["$t"]["load"][l]["qd"][:,n] = load["timeInterval"]["$t"]["timeStep"][string(n-N*(t-1))]["qd"]*sbase
            end
        

            if n > 1 && haskey(resultpf_comp["solution"], "storage")
                for (s,st) in resultpf_comp["solution"]["storage"]
                    case_math["$t"]["storage"][s]["energy"] = st["se"]
                end 
            end

            total_pv = length(collect(keys(optOnlineData["controlled_DER"]["$t"])))
            pv_ind = 0
            pv_ph_ind = 0
            for (k,pv) in optOnlineData["controlled_DER"]["$t"]
                pv_ind += 1
                num_var_pv = 2*length(pv["connections"])
                for (p,ph) in enumerate(pv["connections"])
                    pv_ph_ind += 1 

                    probeSig_p = optOnlineData["probeSig"]["$t"][pv_ph_ind, n+1]
                    probeSig_q = optOnlineData["probeSig"]["$t"][pv_ph_ind+length(pv["connections"]), n+1]
                    explorationResults = OnlineOptDynaGrid.ComputeCostsOPFMATH(optOnlineData, p, ph, t, n-N*(t-1), k, probeSig_p, probeSig_q)
                    @show explorationResults["status_up"], explorationResults["status_down"], k

                    L0_rt_up = Dict{String, Any}()
                    L0_rt_down = Dict{String, Any}()
                    LQ0_rt_up = Dict{String, Any}()
                    LQ0_rt_down = Dict{String, Any}()
                    for (type,data) in optOnlineData["slack"]["$t"]
                        for (g,gen) in data
                            g0_rt_up_slack = vec(explorationResults["g0_rt_up"].*optOnlineData["lam_id"]["$t"][type][g])
                            g0_rt_down_slack = vec(explorationResults["g0_rt_down"].*optOnlineData["lam_id"]["$t"][type][g])

                            L0_rt_up[g] = explorationResults["f0_rt_up"][g] + optOnlineData["lambda_vol_up"]["$t"][g][(ph-1)*2*num_vol["$t"]+1:ph*2*num_vol["$t"],n]'*g0_rt_up_slack
                            L0_rt_down[g] = explorationResults["f0_rt_down"][g] + optOnlineData["lambda_vol_up"]["$t"][g][(ph-1)*2*num_vol["$t"]+1:ph*2*num_vol["$t"],n]'*g0_rt_down_slack
                            LQ0_rt_up[g] = explorationResults["h0_rt_up"][g] + optOnlineData["lambda_vol_up"]["$t"][g][(ph-1)*2*num_vol["$t"]+1:ph*2*num_vol["$t"],n]'*g0_rt_up_slack
                            LQ0_rt_down[g] = explorationResults["h0_rt_down"][g] + optOnlineData["lambda_vol_up"]["$t"][g][(ph-1)*2*num_vol["$t"]+1:ph*2*num_vol["$t"],n]'*g0_rt_down_slack

                            optOnlineData["grad"]["$t"][type][g][pv_ph_ind, n+1]=probeSig_p*(L0_rt_up[g]-L0_rt_down[g])/(2*epsilon)
                            optOnlineData["grad"]["$t"][type][g][pv_ph_ind, n+1]=max.(-0.025, min.(0.025, optOnlineData["grad"]["$t"][type][g][pv_ph_ind, n+1]))   
                            optOnlineData["grad"]["$t"][type][g][pv_ph_ind + length(pv["connections"]), n+1]=probeSig_q*(LQ0_rt_up[g]-LQ0_rt_down[g])/(2*epsilon)
                            optOnlineData["grad"]["$t"][type][g][pv_ph_ind + length(pv["connections"]), n+1]=max.(-0.025, min.(0.025, optOnlineData["grad"]["$t"][type][g][pv_ph_ind + length(pv["connections"]), n+1])) 

                            global optOnlineData["lambda_vol_up"]["$t"][g][(ph-1)*2*num_vol["$t"]+1:ph*2*num_vol["$t"], n+1]=optOnlineData["lambda_vol_up"]["$t"][g][(ph-1)*2*num_vol["$t"]+1:ph*2*num_vol["$t"], n]+ 20*(g0_rt_up_slack+g0_rt_down_slack)/2
                            optOnlineData["lambda_vol_up"]["$t"][g][(ph-1)*2*num_vol["$t"]+1:ph*2*num_vol["$t"], n+1]=max.(0, optOnlineData["lambda_vol_up"]["$t"][g][(ph-1)*2*num_vol["$t"]+1:ph*2*num_vol["$t"], n+1])
                        end
                    end 
        
                    # Only one of the grad matrices should have a non-zero value for each element because each slack/grid-forming DER should only respond to 
                    #    PV connected to the same component and each component should only have one grid-forming DER
                    composite_grad_p = 0
                    composite_grad_q = 0
                    for (type,data) in optOnlineData["slack"]["$t"]
                        for (g,gen) in data
                            composite_grad_p += optOnlineData["grad"]["$t"][type][g][pv_ph_ind, n+1] 
                            composite_grad_q += optOnlineData["grad"]["$t"][type][g][pv_ph_ind+ length(pv["connections"]), n+1] 
                        end
                    end

                    optOnlineData["control"]["$t"][pv_ph_ind, n+1] = -alpha*composite_grad_p
                    optOnlineData["control"]["$t"][pv_ph_ind + length(pv["connections"]), n+1] = -alpha*composite_grad_q

                    # Update the controls one at a time
                    optOnlineData["P_control"]["base"][k][p]= optOnlineData["P_control"]["base"][k][p]+optOnlineData["control"]["$t"][pv_ph_ind, n+1]
                    optOnlineData["P_control"]["base"][k][p]=max.(case_math["$t"]["gen"][k]["pmin"][p], min.(case_math["$t"]["gen"][k]["pmax"][p],optOnlineData["P_control"]["base"][k][p]))
                    optOnlineData["Q_control"]["base"][k][p]= optOnlineData["Q_control"]["base"][k][p] + optOnlineData["control"]["$t"][pv_ph_ind + length(pv["connections"]), n+1]
                    optOnlineData["Q_control"]["base"][k][p]=max.(case_math["$t"]["gen"][k]["qmin"][p], min.(case_math["$t"]["gen"][k]["qmax"][p],optOnlineData["Q_control"]["base"][k][p]))
                    @show optOnlineData["P_control"]["base"][k][p]
                end
                pv_ph_ind += length(pv["connections"])
            end 
        
            global case_comp_math=deepcopy(case_math["$t"])
            for j in keys(pm_baseCase.data["load"])
                case_comp_math["load"][j]["pd"]=pm_baseCase.data["load"][j]["timeInterval"]["$t"]["timeStep"][string(n-N*(t-1))]["pd"]
            end
            
            for (m,pv) in optOnlineData["solar"]["$t"]
                if !haskey(optOnlineData["slack"]["$t"], "solar") || !(m in collect(keys(optOnlineData["slack"]["$t"]["solar"])))
                    case_comp_math["gen"][m]["pmin"]=case_math["$t"]["gen"][m]["pg"]
                    case_comp_math["gen"][m]["pmax"]=case_comp_math["gen"][m]["pmin"]

                    case_comp_math["gen"][m]["qmin"]=case_math["$t"]["gen"][m]["qg"]
                    case_comp_math["gen"][m]["qmax"]=case_comp_math["gen"][m]["qmin"]
                end
            end

            for (m,pv) in optOnlineData["controlled_DER"]["$t"]
                case_comp_math["gen"][m]["pg"]= optOnlineData["P_control"]["base"][m]
                case_comp_math["gen"][m]["pmin"]=case_comp_math["gen"][m]["pg"]
                case_comp_math["gen"][m]["pmax"]=case_comp_math["gen"][m]["pmin"]

                case_comp_math["gen"][m]["qg"]= optOnlineData["Q_control"]["base"][m]
                case_comp_math["gen"][m]["qmin"]=case_comp_math["gen"][m]["qg"]
                case_comp_math["gen"][m]["qmax"]=case_comp_math["gen"][m]["qmin"]
            end

            global resultpf_comp = PMD.solve_mc_model(case_comp_math, PMD.IVRUPowerModel, optOnlineData["NL_solver"], OnlineOptDynaGrid.build_mc_pf_oltc_dg; solution_processors=[PMD.sol_data_model!], make_si=false)
            @show resultpf_comp["termination_status"]
            solpfmath_comp=resultpf_comp["solution"]
            for (b,bus) in solpfmath_comp["bus"] 
                optOnlineData["tracked_data"]["$t"]["bus"][b]["vm"][:,n] = bus["vm"]*case_math["$t"]["bus"][b]["status"]
            end

            # Tracking all slack buses in the network
            ptrack_comp = Dict{String, Any}()
            qtrack_comp = Dict{String, Any}()
            f0_rt_comp = Dict{String, Any}()
            h0_rt_comp = Dict{String, Any}()
            for (type,data) in optOnlineData["slack"]["$t"]
                for (g,gen) in data
                    if type == "storage"
                        ptrack_comp[g] = resultpf_comp["solution"]["storage"][g]["ps"] 
                        qtrack_comp[g] = resultpf_comp["solution"]["storage"][g]["qs"] 
                        f0_rt_comp[g] = sum((1*case_math["$t"]["storage"][g]["ps"]-1*ptrack_comp[g]).^2)
                        h0_rt_comp[g] = sum((1*case_math["$t"]["storage"][g]["qs"]-1*qtrack_comp[g]).^2)
                        optOnlineData["tracked_data"]["$t"]["P_track"][type][g]["ps"][:,n] = ptrack_comp[g]*sbase
                        optOnlineData["tracked_data"]["$t"]["P_track"][type][g]["qs"][:,n] = qtrack_comp[g]*sbase
                    elseif type == "voltage_source"
                        # Power injection control
                        ptrack_comp[g] = resultpf_comp["solution"]["gen"][g]["pg"] 
                        qtrack_comp[g] = resultpf_comp["solution"]["gen"][g]["qg"] 
                        f0_rt_comp[g] = sum((1*case_math["$t"]["gen"][g]["pg"]-1*ptrack_comp[g]).^2)
                        h0_rt_comp[g] = sum((1*case_math["$t"]["gen"][g]["qg"]-1*qtrack_comp[g]).^2)

                        optOnlineData["tracked_data"]["$t"]["P_track"][type][g]["pg"][:,n] = ptrack_comp[g]*sbase
                        optOnlineData["tracked_data"]["$t"]["P_track"][type][g]["qg"][:,n] = qtrack_comp[g]*sbase
                    elseif type == "solar"
                        # Power injection control
                        ptrack_comp[g] = resultpf_comp["solution"]["gen"][g]["pg"] 
                        qtrack_comp[g] = resultpf_comp["solution"]["gen"][g]["qg"] 
                        f0_rt_comp[g] = sum((1*case_math["$t"]["gen"][g]["pg"]-1*ptrack_comp[g]).^2)
                        h0_rt_comp[g] = sum((1*case_math["$t"]["gen"][g]["qg"]-1*qtrack_comp[g]).^2)

                        optOnlineData["tracked_data"]["$t"]["P_track"][type][g]["pg"][:,n] = ptrack_comp[g]*sbase
                        optOnlineData["tracked_data"]["$t"]["P_track"][type][g]["qg"][:,n] = qtrack_comp[g]*sbase
                    end

                    optOnlineData["tracked_data"]["$t"]["tracking_cost"][type][g][n] = f0_rt_comp[g]
                    @show optOnlineData["tracked_data"]["$t"]["tracking_cost"][type][g][n]     
                end
            end 

            for (s,pv) in solpfmath_comp["gen"]
                if s in collect(keys(optOnlineData["solar"]["$t"]))
                    optOnlineData["tracked_data"]["$t"]["solar"][s]["pg"][:,n] = pv["pg"]*sbase
                    optOnlineData["tracked_data"]["$t"]["solar"][s]["qg"][:,n] = pv["qg"]*sbase
                end

                if s in collect(keys(optOnlineData["non_solar"]["$t"]))
                    optOnlineData["tracked_data"]["$t"]["non_solar"][s]["pg"][:,n] = pv["pg"]*sbase
                    optOnlineData["tracked_data"]["$t"]["non_solar"][s]["qg"][:,n] = pv["qg"]*sbase
                end
            end

            for (l,load) in solpfmath_comp["load"]
                optOnlineData["tracked_data"]["$t"]["load"][l]["pd"][:,n] = load["pd"].*sbase.*case_math["$t"]["load"][l]["status"]
                optOnlineData["tracked_data"]["$t"]["load"][l]["qd"][:,n] = load["qd"].*sbase.*case_math["$t"]["load"][l]["status"]
            end

            if haskey(solpfmath_comp, "storage")
                for (s,st) in solpfmath_comp["storage"]
                    optOnlineData["tracked_data"]["$t"]["storage"][s]["ps"][:,n] = st["ps"]*sbase
                    optOnlineData["tracked_data"]["$t"]["storage"][s]["qs"][:,n] = st["qs"]*sbase
                end
            end

            merge!(optOnlineData["tracked_data"]["$t"], Dict("controlled_DER" => optOnlineData["controlled_DER"]["$t"]))

        end
    end

    optOnlineResults=Dict{String,Any}("tracked_data"=>optOnlineData["tracked_data"], "resultopf"=>optOnlineData["resultopf"], "CC"=>optOnlineData["CC"], "CC_solar"=>optOnlineData["CC_solar"])
    return optOnlineResults, optOnlineData
    
end