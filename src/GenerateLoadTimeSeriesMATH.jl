
function GenerateLoadTimeSeriesMATH(optOnlineData::Dict{String, Any})
    case_math = optOnlineData["case_math"];
    T = optOnlineData["TimePeriods"];
    N = optOnlineData["TimeSteps"];
    basecase_math = optOnlineData["pm_baseCase"].data;
    timeInterval=Dict()
    for k=1:T
        timeInterval[string(k)]=k
    end

    timeStep=Dict()
    for k=1:N
        timeStep[string(k)]=k
    end

    for k in keys(basecase_math["load"])
        basecase_math["load"][k]["timeInterval"]=Dict()
        for t in keys(timeInterval)
            basecase_math["load"][k]["timeInterval"][t]=Dict()
        end
    end

    for k in keys(basecase_math["load"])
        for t in keys(timeInterval)
            basecase_math["load"][k]["timeInterval"][t]["timeStep"]=Dict()
            basecase_math["load"][k]["timeInterval"][t]["pd"]=case_math[t]["load"][k]["pd"]*(1+0.00*(rand(1)[1]))
            basecase_math["load"][k]["timeInterval"][t]["qd"]=case_math[t]["load"][k]["qd"]*(1+0.00*(rand(1)[1]))
        end
    end

    for k in keys(basecase_math["load"])
        for t in keys(timeInterval)
            for j in keys(timeStep)
                basecase_math["load"][k]["timeInterval"][t]["timeStep"][j]=Dict()
            end
        end
    end


    for k in keys(basecase_math["load"])
        for t in keys(timeInterval)
            for j in keys(timeStep)
                basecase_math["load"][k]["timeInterval"][t]["timeStep"][j]["pd"] = case_math[t]["load"][k]["pd"]*(1+0.00*(rand(1)[1]))
                basecase_math["load"][k]["timeInterval"][t]["timeStep"][j]["qd"] = case_math[t]["load"][k]["qd"]*(1+0.00*(rand(1)[1]))
            end
        end
    end

    return optOnlineData
end     

# This version of the function works with a case_math with a key and Dict for each time interval
function GenerateLoadTimeSeriesMATH(optOnlineData::Dict{String, Any})
    basecase_math = optOnlineData["pm_baseCase"].data;
    case_math = optOnlineData["case_math"]
    T = optOnlineData["TimePeriods"];
    N = optOnlineData["TimeSteps"];
    timeInterval=Dict()
    for k=1:T
        timeInterval[string(k)]=k
    end

    timeStep=Dict()
    for k=1:N
        timeStep[string(k)]=k
    end

    for k in keys(basecase_math["load"])
        basecase_math["load"][k]["timeInterval"]=Dict()
        for t in keys(timeInterval)
            basecase_math["load"][k]["timeInterval"][t]=Dict()
        end
    end

    for k in keys(basecase_math["load"])
        for t in keys(timeInterval)
            basecase_math["load"][k]["timeInterval"][t]["timeStep"]=Dict()
            basecase_math["load"][k]["timeInterval"][t]["pd"]=case_math[t]["load"][k]["pd"]*(1+0.00*(rand(1)[1]))
            basecase_math["load"][k]["timeInterval"][t]["qd"]=case_math[t]["load"][k]["qd"]*(1+0.00*(rand(1)[1]))
        end
    end

    for k in keys(basecase_math["load"])
        for t in keys(timeInterval)
            for j in keys(timeStep)
                basecase_math["load"][k]["timeInterval"][t]["timeStep"][j]=Dict()
            end
        end
    end


    for k in keys(basecase_math["load"])
        for t in keys(timeInterval)
            for j in keys(timeStep)
             basecase_math["load"][k]["timeInterval"][t]["timeStep"][j]["pd"] = case_math[t]["load"][k]["pd"]*(1+0.00*(rand(1)[1]))
             basecase_math["load"][k]["timeInterval"][t]["timeStep"][j]["qd"] = case_math[t]["load"][k]["qd"]*(1+0.00*(rand(1)[1]))
            end
        end
    end

    return optOnlineData
end     