
function ComputeCostsOPFMATH(optOnlineData::Dict{String, <:Any}, p::Int64, ph::Int64, t::Int64, n::Int64, u::String, probeSig_p::Float64, probeSig_q::Float64)
    case_math = optOnlineData["case_math"]["$t"]
    basecase_math = optOnlineData["pm_baseCase"].data 
    P_control = optOnlineData["P_control"]
    Q_control = optOnlineData["Q_control"]
    Vmax = optOnlineData["Vmax"]
    Vmin = optOnlineData["Vmin"]
    slack = optOnlineData["slack"]["$t"]

    NL_solver = optOnlineData["NL_solver"]

    gen_bus = case_math["gen"][u]["gen_bus"]
    bus_vm = case_math["bus"]["$gen_bus"]["vmag"]
    bus_va = case_math["bus"]["$gen_bus"]["va"]

    case_up_math=deepcopy(case_math)
    for j in keys(case_math["load"])
        case_up_math["load"][j]["pd"]=basecase_math["load"][j]["timeInterval"][string(t)]["timeStep"][string(n)]["pd"]
        case_up_math["load"][j]["qd"]=basecase_math["load"][j]["timeInterval"][string(t)]["timeStep"][string(n)]["qd"]
    end

    case_up_math["gen"][u]["pg"][p] = min(P_control["base"][u][p]+probeSig_p, case_up_math["gen"][u]["pmax"][p])
    case_up_math["gen"][u]["qg"][p] = min(Q_control["base"][u][p]+probeSig_q, case_up_math["gen"][u]["qmax"][p])

    for (s,pv) in optOnlineData["solar"]["$t"]
        if s in keys(optOnlineData["controlled_DER"]["$t"]) && s != u 
            case_up_math["gen"][s]["pmin"] = P_control["base"][s]
            case_up_math["gen"][s]["pmax"] = P_control["base"][s]
            case_up_math["gen"][s]["qmin"] = Q_control["base"][s]
            case_up_math["gen"][s]["qmax"] = Q_control["base"][s]
        elseif !haskey(slack, "solar") || !(s in collect(keys(slack["solar"])))
            case_up_math["gen"][s]["pmin"] = case_up_math["gen"][s]["pg"]
            case_up_math["gen"][s]["pmax"] = case_up_math["gen"][s]["pg"]
            case_up_math["gen"][s]["qmin"] = case_up_math["gen"][s]["qg"]
            case_up_math["gen"][s]["qmax"] = case_up_math["gen"][s]["qg"]
        end    
    end

    for (s,pv) in optOnlineData["non_solar"]["$t"]
        if s in keys(optOnlineData["controlled_DER"]["$t"]) && s != u 
            case_up_math["gen"][s]["pmin"] = P_control["base"][s]
            case_up_math["gen"][s]["pmax"] = P_control["base"][s]
            case_up_math["gen"][s]["qmin"] = Q_control["base"][s]
            case_up_math["gen"][s]["qmax"] = Q_control["base"][s]
        elseif !haskey(slack, "solar") || !(s in collect(keys(slack["solar"])))
            case_up_math["gen"][s]["pmin"] = case_up_math["gen"][s]["pg"]
            case_up_math["gen"][s]["pmax"] = case_up_math["gen"][s]["pg"]
            case_up_math["gen"][s]["qmin"] = case_up_math["gen"][s]["qg"]
            case_up_math["gen"][s]["qmax"] = case_up_math["gen"][s]["qg"]
        end    
    end
    
    resultpf_up = PMD.solve_mc_model(case_up_math, PMD.IVRUPowerModel, NL_solver, OnlineOptDynaGrid.build_mc_pf_oltc_dg; solution_processors=[PMD.sol_data_model!], make_si=false)
    status_up=resultpf_up["termination_status"]

    # Tracking all slack buses in the network
    ptrack_up = Dict{String, Any}()
    f0_rt_up = Dict{String, Any}()
    qtrack_up = Dict{String, Any}()
    h0_rt_up = Dict{String, Any}()
    for (type,data) in slack
        for (g,gen) in data
            #
            if u in optOnlineData["CC_solar"]["$t"][type][g]
                if type == "storage"
                    ptrack_up[g] = resultpf_up["solution"]["storage"][g]["ps"]
                    f0_rt_up[g] = sum((1*case_math[type][g]["ps"][ph]-1*ptrack_up[g][ph]).^2)     # using ph here assumes slacks are always 3-phase
                    qtrack_up[g] = resultpf_up["solution"]["storage"][g]["qs"]
                    h0_rt_up[g] = sum((1*case_math[type][g]["qs"][ph]-1*qtrack_up[g][ph]).^2)   
                else
                    # Power injection control
                    ptrack_up[g] = resultpf_up["solution"]["gen"][g]["pg"]
                    f0_rt_up[g] = sum((1*case_math["gen"][g]["pg"][ph]-1*ptrack_up[g][ph]).^2)  
                    qtrack_up[g] = resultpf_up["solution"]["gen"][g]["qg"]
                    h0_rt_up[g] = sum((1*case_math["gen"][g]["qg"][ph]-1*qtrack_up[g][ph]).^2)  
                end
            else
                # Power injection control
                ptrack_up[g] = resultpf_up["solution"]["gen"][g]["pg"]
                f0_rt_up[g] = 0  
                qtrack_up[g] = resultpf_up["solution"]["gen"][g]["qg"]
                h0_rt_up[g] = 0
            end 
        end
    end 

    solpfmath_up=resultpf_up["solution"]
    Vmag_up=Vector{Float64}()
    for k in keys(solpfmath_up["bus"]) 
        append!(Vmag_up, solpfmath_up["bus"][k]["vm"])
    end
    Vsize=length(Vmag_up)
    g0_rt_up=3*[Vmag_up-Vmax*ones(Vsize); Vmin*ones(Vsize)-Vmag_up] # Cofficient will need tuning for different networks
    ################################################################################    


    case_down_math=deepcopy(case_math)
    for j in keys(case_math["load"])
        case_down_math["load"][j]["pd"]=case_up_math["load"][j]["pd"]
        case_down_math["load"][j]["qd"]=case_up_math["load"][j]["qd"]
    end

    case_down_math["gen"][u]["pg"][p] = max(P_control["base"][u][p]-probeSig_p, case_down_math["gen"][u]["pmin"][p])
    case_down_math["gen"][u]["qg"][p] = max(Q_control["base"][u][p]-probeSig_q, case_down_math["gen"][u]["qmin"][p])


    for (s,pv) in optOnlineData["solar"]["$t"]
        if s in keys(optOnlineData["controlled_DER"]["$t"]) && s != u 
            case_down_math["gen"][s]["pmin"] = P_control["base"][s]
            case_down_math["gen"][s]["pmax"] = P_control["base"][s]
            case_down_math["gen"][s]["qmin"] = Q_control["base"][s]
            case_down_math["gen"][s]["qmax"] = Q_control["base"][s]
        elseif !haskey(slack, "solar") || !(s in collect(keys(slack["solar"])))
            case_down_math["gen"][s]["pmin"] = case_down_math["gen"][s]["pg"]
            case_down_math["gen"][s]["pmax"] = case_down_math["gen"][s]["pmin"]
            case_down_math["gen"][s]["qmin"] = case_down_math["gen"][s]["qg"]
            case_down_math["gen"][s]["qmax"] = case_down_math["gen"][s]["qmin"]   
        end     
    end 

    for (s,pv) in optOnlineData["non_solar"]["$t"]
        if s in keys(optOnlineData["controlled_DER"]["$t"]) && s != u 
            case_down_math["gen"][s]["pmin"] = P_control["base"][s]
            case_down_math["gen"][s]["pmax"] = P_control["base"][s]
            case_down_math["gen"][s]["qmin"] = Q_control["base"][s]
            case_down_math["gen"][s]["qmax"] = Q_control["base"][s]
        elseif !haskey(slack, "solar") || !(s in collect(keys(slack["solar"])))
            case_down_math["gen"][s]["pmin"] = case_down_math["gen"][s]["pg"]
            case_down_math["gen"][s]["pmax"] = case_down_math["gen"][s]["pmin"]
            case_down_math["gen"][s]["qmin"] = case_down_math["gen"][s]["qg"]
            case_down_math["gen"][s]["qmax"] = case_down_math["gen"][s]["qmin"]   
        end     
    end 

    resultpf_down = PMD.solve_mc_model(case_down_math, PMD.IVRUPowerModel, NL_solver, OnlineOptDynaGrid.build_mc_pf_oltc_dg; solution_processors=[PMD.sol_data_model!], make_si=false)
    status_down = resultpf_down["termination_status"]
    
    # Tracking all slack buses in the network
    ptrack_down = Dict{String, Any}()
    f0_rt_down = Dict{String, Any}()
    qtrack_down = Dict{String, Any}()
    h0_rt_down = Dict{String, Any}()
    for (type,data) in slack
        for (g,gen) in data
            if u in optOnlineData["CC_solar"]["$t"][type][g]
                if type == "storage"
                    ptrack_down[g] = resultpf_down["solution"]["storage"][g]["ps"]
                    f0_rt_down[g] = sum((1*case_math["storage"][g]["ps"][ph]-1*ptrack_down[g][ph]).^2)
                    qtrack_down[g] = resultpf_down["solution"]["storage"][g]["qs"]
                    h0_rt_down[g] = sum((1*case_math["storage"][g]["qs"][ph]-1*qtrack_down[g][ph]).^2)     
                else
                    # Power injection control
                    ptrack_down[g] = resultpf_down["solution"]["gen"][g]["pg"] 
                    f0_rt_down[g] = sum((1*case_math["gen"][g]["pg"][ph]-1*ptrack_down[g][ph]).^2) 
                    qtrack_down[g] = resultpf_down["solution"]["gen"][g]["qg"] 
                    h0_rt_down[g] = sum((1*case_math["gen"][g]["qg"][ph]-1*qtrack_down[g][ph]).^2) 
                end 
            else
                # Power injection control
                ptrack_down[g] = resultpf_down["solution"]["gen"][g]["pg"]
                f0_rt_down[g] = 0  
                qtrack_down[g] = resultpf_down["solution"]["gen"][g]["qg"]
                h0_rt_down[g] = 0
            end 
        end 
    end 
    solpfmath_down=resultpf_down["solution"]
    
    Vmag_down=Vector{Float64}()
    for k in keys(solpfmath_down["bus"]) 
        append!(Vmag_down, solpfmath_down["bus"][k]["vm"])
    end
    Vsize=length(Vmag_down)
    g0_rt_down=3*[Vmag_down-Vmax*ones(Vsize); Vmin*ones(Vsize)-Vmag_down] # Cofficient will need tuning for different networks

    explorationResults = Dict{String, Any}("f0_rt_up"=>f0_rt_up, "f0_rt_down"=>f0_rt_down, "h0_rt_up"=>h0_rt_up, "h0_rt_down"=>h0_rt_down,
                 "g0_rt_up"=>g0_rt_up, "g0_rt_down"=>g0_rt_down, "status_up"=>status_up, "status_down"=>status_down)
    return explorationResults
end