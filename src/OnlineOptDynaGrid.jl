module OnlineOptDynaGrid
    using Base: Bool


    #Standard Libs
    import Dates
    import LinearAlgebra
    import SHA
    import Statistics
    import Graphs

    #InfrastructureModels
    import InfrastructureModels as IM
    import InfrastructureModels: ismultinetwork

    #PowerModelsDistribution
    import PowerModelsDistribution as PMD
    import PowerModelsDistribution: ref, var, con, sol, ids, nw_ids, nw_id_default, nws
    import PowerModelsDistribution: AbstractUnbalancedPowerModel, ACRUPowerModel, ACPUPowerModel
    import PowerModelsDistribution: IVRUPowerModel, LPUBFDiagModel, LinDist3FlowPowerModel, NFAUPowerModel
    import PowerModelsDistribution: FOTRUPowerModel, FOTPUPowerModel

    #PowerModelsONM
    import PowerModelsONM as ONM

    #Optimization
    import JuMP
    import JuMP: optimizer_with_attributes
    import Ipopt
    import HiGHS
    
    #files
    import JSON
    import JSONSchema
    import CSV

    #combinations
    #import Combinatorics: combinations

    #logging tools
    import Logging
    #import LoggingExtras

    #Revise
    import Revise

    #other tools
    #import Requires: @require

    #=
    function __init()
        global _LOGGER = Logging.ConsoleLogger(;meta_formatter=PMD._pmd_metafmt)
        global _DEFAULT_LOGGER = Logging.current_logger()

        Logging.global_logger(_LOGGER)
        #OnlineOptDynaGrid.set_log_level!(:Info)

        @require Gurobi="2e9cd046-0924-5485-92f1-d5272153d98b" begin
            global GRB_ENV = Gurobi.Env()
        end

        @require KNITRO="67920dd8-b58e-52a8-8622-53c4cffbe346" begin
            global KN_LMC = KNITRO.LMcontext()
        end
    end
    =#

    include("utils.jl")

    include("prob.jl")
    include("dg_variables.jl")
    include("dg_constraint_template.jl")
    include("dg_constraint.jl")
    include("dg_objective.jl")

    include("GenerateLoadTimeSeriesMATH.jl")

    include("ComputeCostsOPFMATH.jl")

    include("runOptOnline.jl")


end # module
