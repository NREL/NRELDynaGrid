

using PowerModelsDistribution

function initializeData(optOnlineData::Dict{String,Any})
    case_math = optOnlineData["case_math"]
    N = optOnlineData["TimeSteps"]
    T = optOnlineData["TimePeriods"]
    N_sim = optOnlineData["N_sim"]

    optOnlineData["P_control"]=Dict("up"=>Dict(), "down"=>Dict(),"base"=>Dict())
    optOnlineData["Q_control"]=Dict("up"=>Dict(), "down"=>Dict(),"base"=>Dict())

    optOnlineData["num_gridForming"] = Dict{String, Int}("$t"=>0 for t in 1:T);
    optOnlineData["grad"] = Dict{String,Dict{String, Dict{String, Any}}}("$t"=>Dict{String, Dict{String, Any}}() for t in 1:T)
    optOnlineData["slack"] = Dict{String,Dict{String, Dict{String, Any}}}("$t"=>Dict{String, Dict{String, Any}}() for t in 1:T)
    optOnlineData["solar"] = Dict{String, Dict{String, Any}}("$t"=>deepcopy(case_math["$t"]["gen"]) for t in 1:T)
    optOnlineData["non_solar"] = Dict{String, Dict{String, Any}}("$t"=>deepcopy(case_math["$t"]["gen"]) for t in 1:T)
    optOnlineData["num_vol"] = Dict{String, Int}("$t"=>0 for t in 1:T);
    optOnlineData["lambda_vol_up"] = Dict{String, Dict{String, Any}}("$t"=>Dict{String, Any}() for t in 1:T);
    optOnlineData["lam_id"] = Dict{String,Dict{String, Dict{String, Any}}}("$t"=>Dict{String, Any}() for t in 1:T)
    optOnlineData["control"] = Dict{String, Matrix}();
    optOnlineData["omega"] = Dict{String, Matrix}();
    optOnlineData["probeSig"] = Dict{String, Matrix}();
    optOnlineData["tracked_data"] = Dict{String, Dict{String, Any}}("$t"=>Dict{String, Any}() for t in 1:T)
    for t in 1:T
        
        #optOnlineData["solar"]["$t"] = deepcopy(case_math["$t"]["gen"])
        for (i,gen) in case_math["$t"]["gen"]
            type = split(gen["source_id"], ".")[1]
            if type != "solar" || (gen["gen_status"] == PMD.DISABLED || gen["gen_status"] == 0)
                delete!(optOnlineData["solar"]["$t"], i)
            end
        end

        #optOnlineData["non_solar"]["$t"] = deepcopy(case_math["$t"]["gen"])
        for (i,gen) in case_math["$t"]["gen"]
            type = split(gen["source_id"], ".")[1]
            if type != "generator" || (gen["gen_status"] == PMD.DISABLED || gen["gen_status"] == 0)
                delete!(optOnlineData["non_solar"]["$t"], i)
            end
        end

        #optOnlineData["num_vol"]["$t"] = 0
        for (b,bus) in case_math["$t"]["bus"]
            if bus["status"] == PMD.ENABLED || bus["status"] == 1
                optOnlineData["num_vol"]["$t"] += length(bus["terminals"])
            end
        end
        nv=optOnlineData["num_vol"]["$t"]

        #optOnlineData["num_gridForming"]["$t"] = 0
        #global optOnlineData["grad"]["$t"] = Dict{String, Dict{String, Any}}()
        #slack["$t"] = Dict{String, Dict{String, Any}}()
        #lambda_vol_up["$t"] = Dict{String, Any}()
        #lam_id["$t"] = Dict{String, Any}()
        for (i,gen) in case_math["$t"]["gen"] # Note that this will produce an error if a generator that is disabled in T=1 is enabled in a later T
            type = split(gen["source_id"], ".")[1]
            if type == "voltage_source" && (gen["gen_status"] == PMD.ENABLED || gen["gen_status"] == 1)
                optOnlineData["num_gridForming"]["$t"] += 1;
                merge!(optOnlineData["grad"]["$t"], Dict("voltage_source" => Dict(i => Vector{Float64}())));
                merge!(optOnlineData["slack"]["$t"], Dict("voltage_source" => Dict(i => gen))) 
                merge!(optOnlineData["lambda_vol_up"]["$t"], Dict(i => zeros(nv*2*3,N_sim+1)))
                merge!(optOnlineData["lam_id"]["$t"], Dict("voltage_source" => Dict(i => Vector{Int64}())))
            elseif type == "solar" && (gen["gen_status"] == PMD.ENABLED || gen["gen_status"] == 1) && (gen["inverter"] == ONM.GRID_FORMING || gen["inverter"] == 1)
                optOnlineData["num_gridForming"]["$t"] += 1;
                if !haskey(optOnlineData["grad"]["$t"], "solar")
                    merge!(optOnlineData["grad"]["$t"], Dict("solar" => Dict(i => Vector{Float64}())));
                    merge!(optOnlineData["slack"]["$t"], Dict("solar" => Dict(i => gen))) 
                    merge!(optOnlineData["lam_id"]["$t"], Dict("solar" =>Dict(i => Vector{Int64}())))
                else
                    merge!(optOnlineData["grad"]["$t"]["solar"], Dict(i => Vector{Float64}()));
                    merge!(optOnlineData["slack"]["$t"]["solar"], Dict(i => gen))
                    merge!(optOnlineData["lam_id"]["$t"]["solar"], Dict(i => Vector{Int64}()));
                end
                merge!(optOnlineData["lambda_vol_up"]["$t"], Dict(i => zeros(nv*2*3,N_sim+1)))
                gen["pmax"] = Inf.*gen["pmax"]
            elseif type == "generator" && (gen["gen_status"] == PMD.ENABLED || gen["gen_status"] == 1) && (gen["inverter"] == ONM.GRID_FORMING || gen["inverter"] == 1)
                optOnlineData["num_gridForming"]["$t"] += 1;
                if !haskey(optOnlineData["grad"]["$t"], "non_solar")
                    merge!(optOnlineData["grad"]["$t"], Dict("non_solar" => Dict(i => Vector{Float64}())));
                    merge!(optOnlineData["slack"]["$t"], Dict("non_solar" => Dict(i => gen))) 
                    merge!(optOnlineData["lam_id"]["$t"], Dict("non_solar" =>Dict(i => Vector{Int64}())))
                else
                    merge!(optOnlineData["grad"]["$t"]["non_solar"], Dict(i => Vector{Float64}()));
                    merge!(optOnlineData["slack"]["$t"]["non_solar"], Dict(i => gen))
                    merge!(optOnlineData["lam_id"]["$t"]["non_solar"], Dict(i => Vector{Int64}()));
                end
                merge!(optOnlineData["lambda_vol_up"]["$t"], Dict(i => zeros(nv*2*3,N_sim+1)))
                gen["pmax"] = Inf.*gen["pmax"]
            end 
        end
        for (i,st) in case_math["$t"]["storage"] # Note that this will produce an error if a generator that is disabled in T=1 is enabled in a later T
            if (st["status"] == PMD.ENABLED || st["status"] == 1) && (st["inverter"] == ONM.GRID_FORMING || st["inverter"] == 1)
                optOnlineData["num_gridForming"]["$t"] += 1;
                if !haskey(optOnlineData["grad"]["$t"], "storage")
                    merge!(optOnlineData["grad"]["$t"], Dict("storage" => Dict(i => Vector{Float64}())));
                    merge!(optOnlineData["slack"]["$t"], Dict("storage" => Dict(i => st))) 
                    merge!(optOnlineData["lam_id"]["$t"], Dict("storage" =>Dict(i => Vector{Int64}())))
                else
                    merge!(optOnlineData["grad"]["$t"]["storage"], Dict(i => Vector{Float64}()));
                    merge!(optOnlineData["slack"]["$t"]["storage"], Dict(i => st))
                    merge!(optOnlineData["lam_id"]["$t"]["storage"], Dict(i => Vector{Int64}()));
                end
                merge!(optOnlineData["lambda_vol_up"]["$t"], Dict(i => zeros(nv*2*3,N_sim+1)))
            end 
        end

        optOnlineData["tracked_data"]["$t"] = OnlineOptDynaGrid.initializeTrackedValues(optOnlineData, t)
    end 

    return optOnlineData
end 

function initializeTrackedValues(optOnlineData::Dict{String, Any}, t::Int)
    tracked_data = optOnlineData["tracked_data"]["$t"]
    N = optOnlineData["TimeSteps"]
    T = optOnlineData["TimePeriods"]
    case = optOnlineData["case_math"]["$t"]

    if case["data_model"] == PMD.ENGINEERING
        buses = Dict{String, Any}()
        for (b,bus) in case["bus"]
            if 4 in bus["connections"]
                ncnds = length(bus["connections"])-1
            else
                ncnds = length(bus["connections"])
            end 
            merge!(buses, Dict(b => Dict("vm" => zeros(ncnds, N*T))))
        end
        merge!(tracked_data, Dict("bus" => buses))

        switches = Dict{String, Any}()
        for (s,sw) in case["switch"]
            if 4 in sw["f_connections"]
                ncnds = length(sw["f_connections"])-1
            else
                ncnds = length(sw["f_connections"])
            end 
            merge!(switches, Dict(s => Dict("pf" => zeros(ncnds, N*T), "qf" => zeros(ncnds, N*T))))
        end
        merge!(tracked_data, Dict("switch" => switches))

        pvs = Dict{String, Any}()
        for (s,pv) in case["solar"]
            if 4 in pv["connections"]
                ncnds = length(pv["connections"])-1
            else
                ncnds = length(pv["connections"])
            end 
            merge!(pvs, Dict(s => Dict("pg" => zeros(ncnds, N*T), "qg" => zeros(ncnds, N*T))))
        end
        merge!(tracked_data, Dict("solar" => pvs))

        non_pvs = Dict{String, Any}()
        for (s,g) in case["generator"]
            if 4 in g["connections"]
                ncnds = length(g["connections"])-1
            else
                ncnds = length(g["connections"])
            end 
            merge!(non_pvs, Dict(s => Dict("pg" => zeros(ncnds, N*T), "qg" => zeros(ncnds, N*T))))
        end
        merge!(tracked_data, Dict("non_solar" => non_pvs))

        loads = Dict{String, Any}()
        for (l,load) in case["load"]
            if 4 in load["connections"]
                ncnds = length(load["connections"])-1
            else
                ncnds = length(load["connections"])
            end 
            merge!(loads, Dict(l => Dict("pd" => zeros(ncnds, N*T), "qd" => zeros(ncnds, N*T))))
        end
        merge!(tracked_data, Dict("load" => loads))

        storage = Dict{String, Any}()
        for (s,st) in case["storage"]
            if 4 in st["connections"]
                ncnds = length(st["connections"])-1
            else
                ncnds = length(st["connections"])
            end 
            merge!(storage, Dict(s => Dict("ps" => zeros(ncnds, N*T), "qs" => zeros(ncnds, N*T))))
        end
        merge!(tracked_data, Dict("storage" => storage))

    elseif case["data_model"] == PMD.MATHEMATICAL
        buses = Dict{String, Any}()
        for (b,bus) in case["bus"]
            ncnds = length(bus["terminals"])
            merge!(buses, Dict(b => Dict("vm" => zeros(ncnds, N*T))))
        end
        merge!(tracked_data, Dict("bus" => buses))

        switches = Dict{String, Any}()
        for (s,sw) in case["switch"]
            if 4 in sw["f_connections"]
                ncnds = length(sw["f_connections"])-1
            else
                ncnds = length(sw["f_connections"])
            end 
            merge!(switches, Dict(s => Dict("pf" => zeros(ncnds, N*T), "qf" => zeros(ncnds, N*T))))
        end
        merge!(tracked_data, Dict("switch" => switches))

        pvs = Dict{String, Any}()
        for (s,pv) in case["gen"]
            type = split(pv["source_id"], ".")[1]
            if type == "solar"
                if 4 in pv["connections"]
                    ncnds = length(pv["connections"])-1
                else
                    ncnds = length(pv["connections"])
                end 
                merge!(pvs, Dict(s => Dict("pg" => zeros(ncnds, N*T), "qg" => zeros(ncnds, N*T))))
            end
        end
        merge!(tracked_data, Dict("solar" => pvs))

        non_pvs = Dict{String, Any}()
        for (s,g) in case["gen"]
            type = split(g["source_id"], ".")[1]
            if type == "generator"
                if 4 in g["connections"]
                    ncnds = length(g["connections"])-1
                else
                    ncnds = length(g["connections"])
                end 
                merge!(non_pvs, Dict(s => Dict("pg" => zeros(ncnds, N*T), "qg" => zeros(ncnds, N*T))))
            end
        end
        merge!(tracked_data, Dict("non_solar" => non_pvs))

        storage = Dict{String, Any}()
        for (s,st) in case["storage"]
            if 4 in st["connections"]
                ncnds = length(st["connections"])-1
            else
                ncnds = length(st["connections"])
            end 
            merge!(storage, Dict(s => Dict("ps" => zeros(ncnds, N*T), "qs" => zeros(ncnds, N*T))))
        end
        merge!(tracked_data, Dict("storage" => storage))
        
        P_track = Dict{String, Dict{String, Any}}()
        grid_forming_DERs = Dict{String, Dict{String, Any}}()
        for (type,data) in optOnlineData["slack"]["$t"]
            merge!(P_track, Dict(type => Dict{String, Any}()))
            merge!(grid_forming_DERs, Dict(type => Dict{String, Any}()))
            for (g,gen) in data
                if 4 in gen["connections"]
                    ncnds = length(gen["connections"])-1
                else
                    ncnds = length(gen["connections"])
                end 
                if type == "storage"
                    merge!(P_track[type], Dict(g => Dict("ps" => zeros(ncnds, N*T), "qs" => zeros(ncnds, N*T))))
                elseif type == "voltage_source" || type == "solar" || type == "generator"
                    merge!(P_track[type], Dict(g => Dict("pg" => zeros(ncnds, N*T), "qg" => zeros(ncnds, N*T))))
                end 
                merge!(grid_forming_DERs[type], Dict(g => zeros(N*T)))
            end 
        end
        merge!(tracked_data, Dict("P_track" => P_track))
        merge!(tracked_data, Dict("tracking_cost" => grid_forming_DERs))

        loads = Dict{String, Any}()
        for (l,load) in case["load"]
            if 4 in load["connections"]
                ncnds = length(load["connections"])-1
            else
                ncnds = length(load["connections"])
            end 
            merge!(loads, Dict(l => Dict("pd" => zeros(ncnds, N*T), "qd" => zeros(ncnds, N*T))))
        end
        merge!(tracked_data, Dict("load" => loads))

    end 

    return tracked_data
end

function fixDataToSolution(optOnlineData::Dict{String, Any})
    case_math = optOnlineData["case_math"]
    case_math_bounds = optOnlineData["case_math_bounds"]
    mld_results = optOnlineData["results"]
    basecase_math = optOnlineData["pm_baseCase"].data
    T = optOnlineData["TimePeriods"]
    Vmax = optOnlineData["Vmax"]
    Vmin = optOnlineData["Vmin"]
    
    for t in 1:T
        for k in keys(case_math["$t"]["load"])
            case_math_bounds["$t"]["load"][k]["pd"]=basecase_math["load"][k]["timeInterval"]["$t"]["timeStep"]["1"]["pd"]
            case_math_bounds["$t"]["load"][k]["qd"]=basecase_math["load"][k]["timeInterval"]["$t"]["timeStep"]["1"]["qd"]
            case_math_bounds["$t"]["load"][k]["status"] = mld_results["$t"]["solution"]["load"][k]["status"] == DISABLED ? 0 : 1
            case_math["$t"]["load"][k]["status"] = mld_results["$t"]["solution"]["load"][k]["status"] == DISABLED ? 0 : 1
            case_math["$t"]["load"][k]["pd"]=basecase_math["load"][k]["timeInterval"]["$t"]["timeStep"]["1"]["pd"]
            case_math["$t"]["load"][k]["qd"]=basecase_math["load"][k]["timeInterval"]["$t"]["timeStep"]["1"]["qd"]
        end

        for (s,sw) in case_math["$t"]["switch"]
            if s in keys(mld_results["$t"]["solution"]["switch"])
                sw["status"] = mld_results["$t"]["solution"]["switch"][s]["state"] == 0 ? 0 : 1
                case_math_bounds["$t"]["switch"][s]["state"] = mld_results["$t"]["solution"]["switch"][s]["state"] == 0 ? 0 : 1
            end
        end 

        for (s,pv) in case_math["$t"]["gen"]
            type = split(pv["source_id"], ".")[1]
            if (type == "solar" || type == "generator") && (pv["gen_status"] == 1 || pv["gen_status"] == ONM.ENABLED) # ignore disabled
                # Set solar status from solution
                pv["gen_status"] = mld_results["$t"]["solution"]["gen"][s]["status"] == DISABLED ? 0 : 1
                case_math_bounds["$t"]["gen"][s]["status"] = mld_results["$t"]["solution"]["gen"][s]["status"] == DISABLED ? 0 : 1
            
                # Set bus type from solution
                pv_bus = pv["gen_bus"]
                if  mld_results["$t"]["solution"]["gen"][s]["inverter"] == 1 || mld_results["$t"]["solution"]["gen"][s]["inverter"] == ONM.GRID_FORMING
                    case_math["$t"]["bus"]["$pv_bus"]["bus_type"] = 3
                    case_math_bounds["$t"]["bus"]["$pv_bus"]["bus_type"] = 3
                    pv["inverter"] =  1
                    case_math_bounds["$t"]["gen"][s]["inverter"] = 1
                    case_math["$t"]["bus"]["$pv_bus"]["vm"] = [1.0, 1.0, 1.0]
                    case_math_bounds["$t"]["bus"]["$pv_bus"]["vm"] = [1.0, 1.0, 1.0]
                    case_math["$t"]["bus"]["$pv_bus"]["vmin"] = [1.0, 1.0, 1.0]
                    case_math_bounds["$t"]["bus"]["$pv_bus"]["vmin"] = [1.0, 1.0, 1.0]
                    case_math["$t"]["bus"]["$pv_bus"]["vmax"] = [1.0, 1.0, 1.0]
                    case_math_bounds["$t"]["bus"]["$pv_bus"]["vmax"] = [1.0, 1.0, 1.0]
                elseif mld_results["$t"]["solution"]["gen"][s]["inverter"] == 0 || mld_results["$t"]["solution"]["gen"][s]["inverter"] == ONM.GRID_FOLLOWING 
                    pv["inverter"] = 0
                    case_math_bounds["$t"]["gen"][s]["inverter"] = 0
                    case_math["$t"]["bus"]["$pv_bus"]["bus_type"] = 2
                    case_math_bounds["$t"]["bus"]["$pv_bus"]["bus_type"] = 2
                    case_math["$t"]["bus"]["$pv_bus"]["vmin"] = [Vmin-0.1, Vmin-0.1, Vmin-0.1]
                    case_math_bounds["$t"]["bus"]["$pv_bus"]["vmin"] = [Vmin, Vmin, Vmin]
                    case_math["$t"]["bus"]["$pv_bus"]["vmax"] = [Vmax+0.1, Vmax+0.1, Vmax+0.1]
                    case_math_bounds["$t"]["bus"]["$pv_bus"]["vmax"] = [Vmax, Vmax, Vmax]
                    if haskey(case_math["$t"]["bus"]["$pv_bus"], "vm")
                        delete!(case_math["$t"]["bus"]["$pv_bus"], "vm")
                        delete!(case_math["$t"]["bus"]["$pv_bus"], "va")
                    end
                    if haskey(case_math_bounds["$t"]["bus"]["$pv_bus"], "vm")
                        delete!(case_math_bounds["$t"]["bus"]["$pv_bus"], "vm")
                        delete!(case_math_bounds["$t"]["bus"]["$pv_bus"], "va")
                    end
                end  
            end
        end 

        for (s,st) in case_math["$t"]["storage"]
            # Set storage status from solution
            if haskey(mld_results["$t"]["solution"], "storage") && st["status"] == 1
                st["status"] = mld_results["$t"]["solution"]["storage"][s]["status"] == DISABLED ? 0 : 1
                case_math_bounds["$t"]["storage"][s]["status"] = mld_results["$t"]["solution"]["storage"][s]["status"] == DISABLED ? 0 : 1
                st_bus = st["storage_bus"]

                if  mld_results["$t"]["solution"]["storage"][s]["inverter"] == 1 || mld_results["$t"]["solution"]["storage"][s]["inverter"] == ONM.GRID_FORMING
                    st["inverter"] =  1
                    case_math_bounds["$t"]["storage"][s]["inverter"] = 1
                    case_math["$t"]["bus"]["$st_bus"]["bus_type"] = 3
                    case_math_bounds["$t"]["bus"]["$st_bus"]["bus_type"] = 3
                    case_math["$t"]["bus"]["$st_bus"]["vm"] = [1.0, 1.0, 1.0]
                    case_math_bounds["$t"]["bus"]["$st_bus"]["vm"] = [1.0, 1.0, 1.0]
                    case_math["$t"]["bus"]["$st_bus"]["vmin"] = [1.0, 1.0, 1.0]
                    case_math_bounds["$t"]["bus"]["$st_bus"]["vmin"] = [1.0, 1.0, 1.0]
                    case_math["$t"]["bus"]["$st_bus"]["vmax"] = [1.0, 1.0, 1.0]
                    case_math_bounds["$t"]["bus"]["$st_bus"]["vmax"] = [1.0, 1.0, 1.0]
                elseif mld_results["$t"]["solution"]["storage"][s]["inverter"] == 0 || mld_results["$t"]["solution"]["storage"][s]["inverter"] == ONM.GRID_FOLLOWING 
                    st["inverter"] = 0
                    case_math_bounds["$t"]["storage"][s]["inverter"] = 0
                    case_math["$t"]["bus"]["$st_bus"]["bus_type"] = 2
                    case_math_bounds["$t"]["bus"]["$st_bus"]["bus_type"] = 2
                    case_math["$t"]["bus"]["$st_bus"]["vmin"] = [Vmin-0.1, Vmin-0.1, Vmin-0.1]
                    case_math_bounds["$t"]["bus"]["$st_bus"]["vmin"] = [min, Vmin, Vmin]
                    case_math["$t"]["bus"]["$st_bus"]["vmax"] = [Vmax+0.1, Vmax+0.1, Vmax+0.1]
                    case_math_bounds["$t"]["bus"]["$st_bus"]["vmax"] = [Vmax, Vmax, Vmax]
                    if haskey(case_math["$t"]["bus"]["$st_bus"], "vm")
                        delete!(case_math["$t"]["bus"]["$st_bus"], "vm")
                        delete!(case_math["$t"]["bus"]["$st_bus"], "va")
                    end
                    if haskey(case_math_bounds["$t"]["bus"]["$st_bus"], "vm")
                        delete!(case_math_bounds["$t"]["bus"]["$st_bus"], "vm")
                        delete!(case_math_bounds["$t"]["bus"]["$st_bus"], "va")
                    end
                end    
            end
        end 

        for (b,bus) in case_math["$t"]["bus"]
            # Set bus status from solution
            bus["status"] = mld_results["$t"]["solution"]["bus"][b]["status"] == DISABLED ? 0 : 1
            case_math_bounds["$t"]["bus"][b]["status"] = mld_results["$t"]["solution"]["bus"][b]["status"] == DISABLED ? 0 : 1

            if bus["status"] == 0
                bus["bus_type"] = 4
            end  
            if case_math_bounds["$t"]["bus"][b]["status"] == 0
                case_math_bounds["$t"]["bus"][b]["bus_type"] = 4
            end 
        end
    end 

    return optOnlineData
end

function fixGen(optOnlineData::Dict{String,Any}, t::Int)
    case_math = optOnlineData["case_math"]["$t"]
    solar = optOnlineData["solar"]["$t"]
    non_solar = optOnlineData["non_solar"]["$t"]
    controlled_DER = optOnlineData["controlled_DER"]["$t"]
    slack = optOnlineData["slack"]["$t"]
    resultopf = optOnlineData["resultopf"]["$t"]
    P_control = optOnlineData["P_control"]
    Q_control = optOnlineData["Q_control"]
    
    for (s,pv) in solar
        # Set solar limits from solution
        case_math["gen"][s]["pg"] = copy(resultopf["solution"]["gen"][s]["pg"])
        case_math["gen"][s]["qg"] = copy(resultopf["solution"]["gen"][s]["qg"])

        if haskey(resultopf["solution"]["gen"][s], "crg")
            case_math["gen"][s]["crg"] = copy(resultopf["solution"]["gen"][s]["crg"])
        end 
        if haskey(resultopf["solution"]["gen"][s], "cig")
            case_math["gen"][s]["cig"] = copy(resultopf["solution"]["gen"][s]["cig"])
        end 

        if s in collect(keys(controlled_DER))
            P_control["base"][s]=copy(resultopf["solution"]["gen"][s]["pg"])
            Q_control["base"][s]=copy(resultopf["solution"]["gen"][s]["qg"])
        end 
    end
    
    for (s,g) in non_solar
        # Set solar limits from solution
        case_math["gen"][s]["pg"] = copy(resultopf["solution"]["gen"][s]["pg"])
        case_math["gen"][s]["qg"] = copy(resultopf["solution"]["gen"][s]["qg"])

        if haskey(resultopf["solution"]["gen"][s], "crg")
            case_math["gen"][s]["crg"] = copy(resultopf["solution"]["gen"][s]["crg"])
        end 
        if haskey(resultopf["solution"]["gen"][s], "cig")
            case_math["gen"][s]["cig"] = copy(resultopf["solution"]["gen"][s]["cig"])
        end 

        if s in collect(keys(optOnlineData["controlled_DER"]))
            P_control["base"][s]=copy(resultopf["solution"]["gen"][s]["pg"])
            Q_control["base"][s]=copy(resultopf["solution"]["gen"][s]["qg"])
        end 
    end

    for (type,data) in slack
        for (g,gen) in data
        # Set ref limits from solution
            if type == "storage"
                case_math["storage"][g]["ps"] = copy(resultopf["solution"]["storage"][g]["ps"])
                case_math["storage"][g]["qs"] = copy(resultopf["solution"]["storage"][g]["qs"])
            else
                case_math["gen"][g]["pg"] = copy(resultopf["solution"]["gen"][g]["pg"])
                case_math["gen"][g]["qg"] = copy(resultopf["solution"]["gen"][g]["qg"])
                if haskey(resultopf["solution"]["gen"][g], "crg")
                    case_math["gen"][g]["crg"] = copy(resultopf["solution"]["gen"][g]["crg"])
                end 
                if haskey(resultopf["solution"]["gen"][g], "cig")
                    case_math["gen"][g]["cig"] = copy(resultopf["solution"]["gen"][g]["cig"])
                end 
            end 
        end
    end

    for (s,st) in case_math["storage"]
        # Set storage status from solution
        if haskey(resultopf["solution"], "storage") && st["status"] == 1
            st["ps"] = resultopf["solution"]["storage"][s]["ps"]
            st["qs"] = resultopf["solution"]["storage"][s]["qs"]
        end
    end 

    for (b,bus) in case_math["bus"]
        if bus["status"] == 1
            bus["vmag"] = copy(resultopf["solution"]["bus"][b]["vm"])
            bus["va"] = copy(resultopf["solution"]["bus"][b]["va"])
            bus["vm_start"] = copy(resultopf["solution"]["bus"][b]["vm"])
            bus["va_start"] = copy(resultopf["solution"]["bus"][b]["va"])
        end 
    end 

    return optOnlineData
end

function fill_lam_id(optOnlineData::Dict{String, Any})
    lam_id = optOnlineData["lam_id"]
    pm = optOnlineData["pm_baseCase"]
    case_math = optOnlineData["case_math"]
    N_sim = optOnlineData["N_sim"]
    grad = optOnlineData["grad"]
    resultopf = optOnlineData["resultopf"]
    
    map_id_pairs = Dict{Int,Tuple{Int,Int}}(id => (PMD.ref(pm, 0, :bus_block_map, sw["f_bus"]),PMD.ref(pm, 0, :bus_block_map, sw["t_bus"])) for (id,sw) in PMD.ref(pm, 0, :switch))
    block_sws = PMD.ref(pm, 0, :block_switches)
    sbase = case_math["1"]["settings"]["sbase_default"]

    # create graph of blocks
    graph = Graphs.SimpleGraph(length(keys(PMD.ref(pm, 0, :blocks))))
    for (s, (b1,b2)) in map_id_pairs
        if b1 != b2 
            Graphs.add_edge!(graph, b1, b2)
        end
    end

    optOnlineData["CC"] = Dict{String, Dict{String, Dict{String, Vector{String}}}}(t => Dict(type => Dict(gfm => Vector{String}() for (gfm,_) in vals) for (type,vals) in data) for (t,data) in lam_id)
    optOnlineData["CC_solar"] = Dict{String, Dict{String, Dict{String, Vector{String}}}}(t => Dict(type => Dict(gfm => Vector{String}() for (gfm,_) in vals) for (type,vals) in data) for (t,data) in lam_id)
    optOnlineData["controlled_DER"] = Dict{String, Dict{String, Any}}(t => Dict{String, Any}() for (t,data) in case_math)
    optOnlineData["N_var"] = Dict{String, Int}(t => 0 for (t,data) in case_math);
    for (t,data) in lam_id
        graph_t = copy(graph)
        for (s,sw) in case_math[t]["switch"]
            if sw["status"] == 0 && haskey(map_id_pairs, parse(Int,s))
                Graphs.rem_edge!(graph_t, map_id_pairs[parse(Int,s)][1], map_id_pairs[parse(Int,s)][2])
            end
        end
        for (type,vals) in data
            for (gfm,_) in vals
                gfm_int = parse(Int, gfm)
                gfm_block = PMD.ref(pm, 0, :gen_block_map)[gfm_int]
                for (bl,buses) in PMD.ref(pm, 0, :blocks)
                    if Graphs.has_path(graph_t, gfm_block, bl)
                        for bus in buses
                            append!(optOnlineData["CC"][t][type][gfm], ["$bus"])
                        end
                    end
                end

                block_phs = case_math[t]["gen"][gfm]["connections"]
                largest_pv_per_ph = Dict{Int, Tuple}(ph => ("", 0.0) for ph in block_phs)
                second_largest = Dict{Int, Tuple}(ph => ("", 0.0) for ph in block_phs)
                for (g,gen) in case_math[t]["gen"]
                    if (string(gen["gen_bus"]) in optOnlineData["CC"][t][type][gfm]) && (g != gfm)
                        if haskey(resultopf[t]["solution"]["gen"], g)
                            append!(optOnlineData["CC_solar"][t][type][gfm], [g])

                            for ph in gen["connections"] 
                                if ph in block_phs
                                    room = gen["pmax"] - resultopf[t]["solution"]["gen"][g]["pg"]
                                    if minimum(room) > largest_pv_per_ph[ph][2]
                                        second_largest[ph] = largest_pv_per_ph[ph]
                                        largest_pv_per_ph[ph] = (g,minimum(room))
                                    end
                                end
                            end
                        end 
                    end
                end 

                for ph in block_phs
                    if second_largest[ph] == ("", 0.0) || second_largest[ph][1] == largest_pv_per_ph[ph][1]
                        second_largest[ph] = ("", 0.0)
                        for (g,gen) in case_math[t]["gen"]
                            if ((string(gen["gen_bus"]) in optOnlineData["CC"][t][type][gfm]) && (g != gfm))  
                                if haskey(resultopf[t]["solution"]["gen"], g) && (g != largest_pv_per_ph[ph][1])
                                    for ph in gen["connections"] 
                                        if ph in block_phs
                                            room = gen["pmax"] - resultopf[t]["solution"]["gen"][g]["pg"]
                                            if minimum(room) > second_largest[ph][2]
                                                second_largest[ph] = (g,minimum(room))
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end 
            
                for (ph, (pv,pmax)) in largest_pv_per_ph
                    merge!(optOnlineData["controlled_DER"][t], Dict(pv => case_math[t]["gen"][pv]))
                    optOnlineData["N_var"][t] += 2*length(case_math[t]["gen"][pv]["connections"])

                    if maximum(case_math[t]["gen"][pv]["pmax"]) < 10/sbase && length(optOnlineData["CC_solar"][t][type][gfm]) > 1
                        second_solar = second_largest[ph]
                        merge!(optOnlineData["controlled_DER"][t], Dict(second_solar[1] => case_math[t]["gen"][second_solar[1]]))
                        optOnlineData["N_var"][t] += 2*length(case_math[t]["gen"][second_solar[1]]["connections"])

                        total_room_up = [minimum(case_math[t]["gen"][pv]["pmax"] - resultopf[t]["solution"]["gen"][pv]["pg"]), minimum(case_math[t]["gen"][second_largest[ph][1]]["pmax"] - resultopf[t]["solution"]["gen"][second_largest[ph][1]]["pg"])]
                        room_up = Dict{String, Float64}()
                        total_room_down = [minimum(resultopf[t]["solution"]["gen"][pv]["pg"] - case_math[t]["gen"][pv]["pmin"]), minimum(resultopf[t]["solution"]["gen"][second_largest[ph][1]]["pg"] - case_math[t]["gen"][second_largest[ph][1]]["pmin"])]
                        room_down = Dict{String, Float64}()
                        for g in optOnlineData["CC_solar"][t][type][gfm] 
                            gen = case_math[t]["gen"][g]
                            if gen["gen_status"] == 1 && (g != gfm)
                                for (ind, phase) in enumerate(gen["connections"])
                                    if phase == ph && g != second_largest[ph][1] && g != largest_pv_per_ph[ph][1]
                                        room_up[g] = gen["pmax"][ind] - resultopf[t]["solution"]["gen"][g]["pg"][ind]
                                        room_down[g] = resultopf[t]["solution"]["gen"][g]["pg"][ind] - gen["pmin"][ind]
                                    end 
                                end
                            end 
                        end

                        while sum(total_room_up) < 10/sbase && !isempty(room_up)
                            (val_up, ind_up) = findmax(room_up)
                            merge!(optOnlineData["controlled_DER"][t], Dict(ind_up => case_math[t]["gen"][ind_up]))
                            optOnlineData["N_var"][t] += 2*length(case_math[t]["gen"][ind_up]["connections"])
                            push!(total_room_up, minimum(case_math[t]["gen"][ind_up]["pmax"] - resultopf[t]["solution"]["gen"][ind_up]["pg"]))
                            delete!(room_up, ind_up)

                        end 
                        while sum(total_room_down) < 10/sbase 
                            (val_down, ind_down) = findmax(room_down)
                            merge!(optOnlineData["controlled_DER"][t], Dict(ind_down => case_math[t]["gen"][ind_down]))
                            optOnlineData["N_var"][t] += 2*length(case_math[t]["gen"][ind_down]["connections"])
                            push!(total_room_down, minimum(resultopf[t]["solution"]["gen"][ind_down]["pg"] - case_math[t]["gen"][ind_down]["pmin"]))
                            delete!(room_down, ind_down)
                        end 
                    end 
                end 
             

                for (b,bus) in case_math[t]["bus"]
                    if b in optOnlineData["CC"][t][type][gfm] && bus["status"] == 1
                        append!(lam_id[t][type][gfm], ones(Int, 2*length(bus["terminals"])))
                    elseif !(b in optOnlineData["CC"][t][type][gfm]) && bus["status"] == 1
                        append!(lam_id[t][type][gfm], zeros(Int, 2*length(bus["terminals"])))
                    end
                end 
            end 
        end 
    end 

    for (t,data) in lam_id
        for (type, vals) in data
            for (gfm,_) in vals
                grad[t][type][gfm] = zeros(Float64, (optOnlineData["N_var"][t], N_sim+1))
            end 
        end
    end 

    return optOnlineData
end

function sampleProbeSignal(optOnlineData)
    N = optOnlineData["TimeSteps"]
    T = optOnlineData["TimePeriods"]
    N_var = optOnlineData["N_var"]
    epsilon = optOnlineData["epsilon"]
    N_sim = optOnlineData["N_sim"]

    for t in 1:T
        optOnlineData["probeSig"]["$t"] = zeros(N_var["$t"], N_sim);
        optOnlineData["omega"]["$t"] = zeros(N_var["$t"], 1)
        for i = 1:N_var["$t"]
            optOnlineData["omega"]["$t"][i]=1. /(optOnlineData["a0"]*i+optOnlineData["a1"]) * 2 * pi; # decreasing the coeff after + in () makes oscillations
            for b = 1:N
                optOnlineData["probeSig"]["$t"][i,(N*t)-(N-b)] = epsilon*(cos(optOnlineData["omega"]["$t"][i].*b));
            end
        end 
        optOnlineData["control"]["$t"] = zeros(N_var["$t"], N_sim+1);
    end 

    return optOnlineData
end