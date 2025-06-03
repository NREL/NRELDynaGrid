
# # # # # # # # # # Constraints for IVR form # # # # # # # # # # # # # #
"""
    constraint_mc_theta_ref(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for reference angle constraints.
"""
function optOn_constraint_mc_theta_ref(pm::PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    bus = PMD.ref(pm, nw, :bus, i)
    terminals = bus["terminals"]
    if haskey(bus, "va")
        va_ref = get(PMD.ref(pm, nw, :bus, i), "va", [deg2rad.([0.0, -120.0, 120.0])..., zeros(length(terminals))...][terminals])
        optOn_constraint_mc_theta_ref(pm, nw, i, va_ref)
    end
    nothing
end

function optOn_constraint_mc_storage_losses(pm::PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    storage = PMD.ref(pm, nw, :storage, i)

    optOn_constraint_mc_storage_losses(pm, nw, i, storage["storage_bus"], storage["connections"], storage["r"], storage["x"], storage["p_loss"], storage["q_loss"])
    nothing
end

function optOn_constraint_mc_storage_losses(pm::PMD.AbstractUnbalancedIVRModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    storage = PMD.ref(pm, nw, :storage, i)

    optOn_constraint_mc_storage_losses(pm, nw, i, storage["storage_bus"], storage["connections"], storage["r"], storage["x"], storage["p_loss"], storage["q_loss"])
    nothing
end

"""
    constraint_storage_state(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage state constraints (non multinetwork)
"""
function optOn_constraint_storage_state(pm::PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    storage = PMD.ref(pm, nw, :storage, i)

    if haskey(PMD.ref(pm, nw), :time_elapsed)
        time_elapsed = PMD.ref(pm, nw, :time_elapsed)
    else
        @warn "network data should specify time_elapsed in hours, using 1.0 as a default"
        time_elapsed = 1.0
    end
    
    optOn_constraint_storage_state_initial(pm, nw, i, storage["energy"], storage["charge_efficiency"], storage["discharge_efficiency"], time_elapsed, storage["connections"])
    nothing
end 

function optOn_constraint_storage_state(pm::PMD.AbstractUnbalancedPowerModel, i::Int, nw_1::Int, nw_2::Int)::Nothing
    storage = PMD.ref(pm, nw_2, :storage, i)

    if haskey(PMD.ref(pm, nw_2), :time_elapsed)
        time_elapsed = PMD.ref(pm, nw_2, :time_elapsed)
    else
        @warn "network $(nw_2) should specify time_elapsed in hours, using 1.0 as a default"
        time_elapsed = 1.0
    end

    if haskey(PMD.ref(pm, nw_1, :storage), i)
        optOn_constraint_storage_state(pm, nw_1, nw_2, i, storage["charge_efficiency"], storage["discharge_efficiency"], time_elapsed, storage["connections"])
    else
        println("ERROR: Missing storage state")
    end
    nothing
end

"""
    constraint_mc_storage_power_setpoint_real(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage active power setpoint constraint, for power flow problems
"""
function constraint_mc_storage_power_setpoint_real(pm::AbstractUnbalancedIVRModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    ps_set = ref(pm, nw, :storage, i)["ps"]
    constraint_mc_storage_power_setpoint_real(pm, nw, i, ps_set)
    nothing
end

"""
    constraint_mc_storage_power_setpoint_real(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage active power setpoint constraint, for power flow problems
"""
function constraint_mc_storage_power_setpoint_imaginary(pm::AbstractUnbalancedIVRModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    qs_set = ref(pm, nw, :storage, i)["qs"]
    constraint_mc_storage_power_setpoint_imaginary(pm, nw, i, qs_set)
    nothing
end

"""
    constraint_mc_storage_power_setpoint_real(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage active power setpoint constraint, for power flow problems
"""
function constraint_mc_gen_power_setpoint_imaginary(pm::AbstractUnbalancedIVRModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    qg_set = ref(pm, nw, :gen, i)["qg"]
    constraint_mc_gen_power_setpoint_imaginary(pm, nw, i, qg_set)
    nothing
end

function constraint_mc_gen_current_setpoint_real(pm::AbstractUnbalancedIVRModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    if haskey(ref(pm, nw, :gen, i),"crg")
        crg_set = ref(pm, nw, :gen, i)["crg"]
        constraint_mc_gen_current_setpoint_real(pm, nw, i, crg_set)
    end
    nothing
end

function constraint_mc_gen_current_setpoint_imaginary(pm::AbstractUnbalancedIVRModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    if haskey(ref(pm, nw, :gen, i),"cig")
        cig_set = ref(pm, nw, :gen, i)["cig"]
        constraint_mc_gen_current_setpoint_imaginary(pm, nw, i, cig_set)
    end 
    nothing
end

# # # # # # # # # # for IVR form with ideal (lossless) storage # # # # # # # # # # # # # #
function optOn_constraint_mc_storage_losses_idealStorage(pm::PMD.AbstractUnbalancedIVRModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    storage = PMD.ref(pm, nw, :storage, i)

    optOn_constraint_mc_storage_losses_idealStorage(pm, nw, i, storage["storage_bus"], storage["connections"], storage["r"], storage["x"], storage["p_loss"], storage["q_loss"])
    nothing
end

"""
    constraint_storage_state(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for storage state constraints (non multinetwork)
"""
function optOn_constraint_storage_state_idealStorage(pm::PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    storage = PMD.ref(pm, nw, :storage, i)

    if haskey(PMD.ref(pm, nw), :time_elapsed)
        time_elapsed = PMD.ref(pm, nw, :time_elapsed)
    else
        @warn "network data should specify time_elapsed in hours, using 1.0 as a default"
        time_elapsed = 1.0
    end
    
    optOn_constraint_storage_state_initial_idealStorage(pm, nw, i, storage["energy"], time_elapsed, storage["connections"])
    nothing
end 

function optOn_constraint_storage_state_idealStorage(pm::PMD.AbstractUnbalancedPowerModel, i::Int, nw_1::Int, nw_2::Int)::Nothing
    storage = PMD.ref(pm, nw_2, :storage, i)

    if haskey(PMD.ref(pm, nw_2), :time_elapsed)
        time_elapsed = PMD.ref(pm, nw_2, :time_elapsed)
    else
        @warn "network $(nw_2) should specify time_elapsed in hours, using 1.0 as a default"
        time_elapsed = 1.0
    end

    if haskey(PMD.ref(pm, nw_1, :storage), i)
        optOn_constraint_storage_state_idealStorage(pm, nw_1, nw_2, i, time_elapsed, storage["connections"])
    else
        println("ERROR: Missing storage state")
    end
    nothing
end