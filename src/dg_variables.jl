""
function variable_mc_enabled_bus_voltage(pm::PMD.AbstractUnbalancedACRModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    variable_mc_enabled_bus_voltage_real(pm; nw=nw, bounded=bounded, report=report)
    variable_mc_enabled_bus_voltage_imaginary(pm; nw=nw, bounded=bounded, report=report)

    # local infeasbility issues without proper initialization;
    # convergence issues start when the equivalent angles of the starting point
    # are further away than 90 degrees from the solution (as given by ACP)
    # this is the default behaviour of _PM, initialize all phases as (1,0)
    # the magnitude seems to have little effect on the convergence (>0.05)
    # updating the starting point to a balanced phasor does the job
    for id in PMD.ids(pm, nw, :bus)
        busref = PMD.ref(pm, nw, :bus, id)
        status = busref["status"]
        if status == 1
            terminals = busref["terminals"]
            grounded = busref["grounded"]

            ncnd = length(terminals)

            if haskey(busref, "vr_start") && haskey(busref, "vi_start")
                vr = busref["vr_start"]
                vi = busref["vi_start"]
            else
                vm_start = fill(1.0, 3)
                for t in 1:3
                    if t in terminals
                        vmax = busref["vmax"][findfirst(isequal(t), terminals)]
                        vm_start[t] = min(vm_start[t], vmax)

                        vmin = busref["vmin"][findfirst(isequal(t), terminals)]
                        vm_start[t] = max(vm_start[t], vmin)
                    end
                end

                vm = haskey(busref, "vm_start") ? busref["vm_start"] : haskey(busref, "vm") ? busref["vm"] : [vm_start..., fill(0.0, ncnd)...][terminals]
                va = haskey(busref, "va_start") ? busref["va_start"] : haskey(busref, "va") ? busref["va"] : [deg2rad.([0, -120, 120])..., zeros(length(terminals))...][terminals]

                vr = vm .* cos.(va)
                vi = vm .* sin.(va)
            end

            for (idx,t) in enumerate(terminals)
                JuMP.set_start_value(PMD.var(pm, nw, :vr, id)[t], vr[idx])
                JuMP.set_start_value(PMD.var(pm, nw, :vi, id)[t], vi[idx])
            end
        end 
    end

    # apply bounds if bounded
    if bounded
        for i in PMD.ids(pm, nw, :bus)
            if PMD.ref(pm, nw, :bus, i)["status"] == 1
                PMD.constraint_mc_voltage_magnitude_bounds(pm, i; nw=nw)
            end 
        end
    end
end

""
function variable_mc_enabled_bus_voltage_real(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i, bus) in PMD.ref(pm, nw, :bus))

    vr = PMD.var(pm, nw)[:vr] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vr_$(i)",
            start = PMD.comp_start_value(PMD.ref(pm, nw, :bus, i), "vr_start", t, 1.0)
        ) for i in PMD.ids(pm, nw, :bus)
    )

    if bounded
        for (i,bus) in PMD.ref(pm, nw, :bus)
            if haskey(bus, "vmax")
                for (idx,t) in enumerate(terminals[i])
                    PMD.set_lower_bound(vr[i][t], -bus["vmax"][idx])
                    PMD.set_upper_bound(vr[i][t],  bus["vmax"][idx])
                end
            end
        end
    end

    report && IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :vr, PMD.ids(pm, nw, :bus), vr)
end

# # # # # # # # # # # # # # Storage # # # # # # # # # # # # # # 
"variables for modeling storage units, includes grid injection and internal variables"
function optOn_variable_mc_storage_power(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    optOn_variable_mc_storage_power_real(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_mc_storage_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_mc_storage_power_control_imaginary(pm; nw=nw, bounded=bounded, report=report)
    PMD.variable_mc_storage_current(pm; nw=nw, bounded=bounded, report=report)
    PMD.variable_storage_energy(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_storage_charge(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_storage_discharge(pm; nw=nw, bounded=bounded, report=report)
end


"Create variables for `active` storage injection"
function optOn_variable_mc_storage_power_real(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in PMD.ref(pm, nw, :storage))
    ps = PMD.var(pm, nw)[:ps] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_ps_$(i)",
        start = PMD.comp_start_value(PMD.ref(pm, nw, :storage, i), "ps_start", c, 0.0)
    ) for i in PMD.ids(pm, nw, :storage))

    if bounded
        inj_lb, inj_ub = optOn_ref_calc_storage_injection_bounds( PMD.ref(pm, nw, :storage),PMD.ref(pm, nw, :bus))
        for (i, strg) in PMD.ref(pm, nw, :storage)
            for (idx, c) in enumerate(connections[i])
                PMD.set_lower_bound(ps[i][c], min(inj_lb[i][idx], 0.0))
                PMD.set_upper_bound(ps[i][c], max(inj_ub[i][idx], 0.0))
            end
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :storage, :ps, PMD.ids(pm, nw, :storage), ps)
end

"""
   PMD.ref_calc_storage_injection_bounds(storage, buses)

Computes storage bounds
"""
function optOn_ref_calc_storage_injection_bounds(storage, buses)

    injection_lb = Dict()
    injection_ub = Dict()

    for (i, strg) in storage
        connections = strg["connections"]
        if 4 in connections
            ncnds = length(connections)-1
        else
            ncnds = length(connections)
        end 
        injection_lb[i] = fill(-Inf, ncnds)
        injection_ub[i] = fill( Inf, ncnds)

        if haskey(strg, "thermal_rating")
            injection_lb[i] = max.(injection_lb[i], -strg["thermal_rating"]/ncnds)
            injection_ub[i] = min.(injection_ub[i],  strg["thermal_rating"]/ncnds)
        end

        if haskey(strg, "current_rating")
            for (j, t) in connections
                vmax = buses[strg["storage_bus"]]["vmax"][findfirst(isequal(t), buses[strg["storage_bus"]]["terminals"])]

                injection_lb[i][j] = max(injection_lb[i][j], -strg["current_rating"][j]*vmax)
                injection_ub[i][j] = min(injection_ub[i][j],  strg["current_rating"][j]*vmax)
            end
        end
    end

    return injection_lb, injection_ub
end


"Create variables for `reactive` storage injection"
function optOn_variable_mc_storage_power_imaginary(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in PMD.ref(pm, nw, :storage))
    qs = PMD.var(pm, nw)[:qs] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_qs_$(i)",
        start = PMD.comp_start_value(PMD.ref(pm, nw, :storage, i), "qs_start", c, 0.0)
    ) for i in PMD.ids(pm, nw, :storage))

    if bounded
        for (i, strg) in PMD.ref(pm, nw, :storage)
            if 4 in connections[i]
                ncnds = length(connections[i])-1
            else
                ncnds = length(connections[i])
            end 
            if haskey(strg, "qmin")
                for (idx, c) in enumerate(connections[i])
                    PMD.set_lower_bound(qs[i][c], min(strg["qmin"]/ncnds, 0.0))
                end
            end

            if haskey(strg, "qmax")
                for (idx, c) in enumerate(connections[i])
                    PMD.set_upper_bound(qs[i][c], max(strg["qmax"]/ncnds, 0.0))
                end
            end
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :storage, :qs, PMD.ids(pm, nw, :storage), qs)
end

"""
a reactive power slack variable that enables the storage device to inject or
consume reactive power at its connecting bus, subject to the injection limits
of the device.
"""
function optOn_variable_mc_storage_power_control_imaginary(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in PMD.ref(pm, nw, :storage))
    qsc = PMD.var(pm, nw)[:qsc] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_qsc_$(i)",
        start = PMD.comp_start_value(PMD.ref(pm, nw, :storage, i), "qsc_start", c)
    ) for i in PMD.ids(pm, nw, :storage))
    
    if bounded
        inj_lb, inj_ub = optOn_ref_calc_storage_injection_bounds( PMD.ref(pm, nw, :storage), PMD.ref(pm, nw, :bus))
        for (i,storage) in PMD.ref(pm, nw, :storage)
            for (idx, c) in enumerate(connections[i])
                if !isinf(sum(inj_lb[i][idx])) || haskey(storage, "qmin")
                    lb = max(inj_lb[i][idx], get(storage, "qmin", -Inf))
                    PMD.set_lower_bound(qsc[i][c], min(lb, 0.0))
                end
                if !isinf(sum(inj_ub[i][idx])) || haskey(storage, "qmax")
                    ub = min(inj_ub[i][idx], get(storage, "qmax", Inf))
                    PMD.set_upper_bound(qsc[i][c], max(ub, 0.0))
                end
            end
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :storage, :qsc, PMD.ids(pm, nw, :storage), qsc)
end

""
function optOn_variable_mc_storage_charge(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in PMD.ref(pm, nw, :storage))
    sc = PMD.var(pm, nw)[:sc] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_sc_$(i)",
        start = PMD.comp_start_value(PMD.ref(pm, nw, :storage, i), ["sc_start", "sc"], c, 0.0)
    ) for i in PMD.ids(pm, nw, :storage))

    for (i, strg) in PMD.ref(pm, nw, :storage)
        if 4 in connections[i]
            ncnds = length(connections[i])-1
        else
            ncnds = length(connections[i])
        end 

        for (idx, c) in enumerate(connections[i])
            PMD.set_lower_bound(sc[i][c], 0.0)
        end
        if bounded
            if haskey(strg, "charge_rating")
                for (idx, c) in enumerate(connections[i])
                    PMD.set_upper_bound(sc[i][c], max(strg["charge_rating"]/ncnds, 0.0))
                end
            end
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :storage, :sc, PMD.ids(pm, nw, :storage), sc)
end


""
function optOn_variable_mc_storage_discharge(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in PMD.ref(pm, nw, :storage))
    sd = PMD.var(pm, nw)[:sd] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_sd_$(i)",
        start = PMD.comp_start_value(PMD.ref(pm, nw, :storage, i), ["sd_start", "sd"], c, 0.0)
    ) for i in PMD.ids(pm, nw, :storage))
    
    for (i, strg) in PMD.ref(pm, nw, :storage)
        if 4 in connections[i]
            ncnds = length(connections[i])-1
        else
            ncnds = length(connections[i])
        end 

        for (idx, c) in enumerate(connections[i])
            PMD.set_lower_bound(sd[i][c], 0.0)
        end

        if bounded
            if haskey(strg, "charge_rating")
                for (idx, c) in enumerate(connections[i])
                    PMD.set_upper_bound(sd[i][c], max(strg["discharge_rating"]/ncnds, 0.0))
                end
            end
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :storage, :sd, PMD.ids(pm, nw, :storage), sd)
end

# # # # # # # # # # Variables for IVR form # # # # # # # # # # # # # #
# Modified version of the variable_mc_transformer_current() function from PMD, removed code that stores expressions in rectangular power variable space
function optOn_variable_mc_transformer_current(pm::PMD.AbstractUnbalancedIVRModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    PMD.variable_mc_transformer_current_real(pm; nw=nw, bounded=bounded, report=report)
    PMD.variable_mc_transformer_current_imaginary(pm; nw=nw, bounded=bounded, report=report)
end


"variables for modeling storage units, includes grid injection and internal variables"
function optOn_variable_mc_storage_power(pm::PMD.AbstractUnbalancedIVRModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    optOn_variable_mc_storage_power_real(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_mc_storage_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_mc_storage_power_control_imaginary(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_mc_storage_current(pm; nw=nw, bounded=bounded, report=report)
    PMD.variable_storage_energy(pm; nw=nw, bounded=bounded, report=report)

    optOn_variable_mc_storage_charge(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_mc_storage_discharge(pm; nw=nw, bounded=bounded, report=report)
end

function optOn_variable_mc_storage_current(pm::PMD.AbstractUnbalancedIVRModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    optOn_variable_mc_storage_current_real(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_mc_storage_current_imaginary(pm; nw=nw, bounded=bounded, report=report)
end

"variable: `crs[j]` for `j` in `storage`"
function optOn_variable_mc_storage_current_real(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => storage["connections"] for (i,storage) in PMD.ref(pm, nw, :storage))
    crs = PMD.var(pm, nw)[:crs] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_crs_$(i)",
            start = PMD.comp_start_value(PMD.ref(pm, nw, :storage, i), "crs_start", c, 0.0)
        ) for i in PMD.ids(pm, nw, :storage)
    )

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :storage, :crs, PMD.ids(pm, nw, :storage), crs)
end

"variable: `cis[j]` for `j` in `storage`"
function optOn_variable_mc_storage_current_imaginary(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => storage["connections"] for (i,storage) in PMD.ref(pm, nw, :storage))
    cis = PMD.var(pm, nw)[:cis] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_cis_$(i)",
            start = PMD.comp_start_value(PMD.ref(pm, nw, :storage, i), "cis_start", c, 0.0)
        ) for i in PMD.ids(pm, nw, :storage)
    )

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :storage, :cis, PMD.ids(pm, nw, :storage), cis)
end

# # # # # # # # # # for IVR form with ideal (lossless) storage # # # # # # # # # # # # # #
"variables for modeling storage units, includes grid injection and internal variables"
function optOn_variable_mc_storage_power_idealStorage(pm::PMD.AbstractUnbalancedIVRModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    optOn_variable_mc_storage_power_real(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_mc_storage_power_imaginary(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_mc_storage_power_control_imaginary(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_mc_storage_current(pm; nw=nw, bounded=bounded, report=report)
    PMD.variable_storage_energy(pm; nw=nw, bounded=bounded, report=report)

    optOn_variable_mc_storage_charge_idealStorage(pm; nw=nw, bounded=bounded, report=report)
end

""
function optOn_variable_mc_storage_charge_idealStorage(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => strg["connections"] for (i,strg) in PMD.ref(pm, nw, :storage))
    sc = PMD.var(pm, nw)[:sc] = Dict(i => JuMP.@variable(pm.model,
        [c in connections[i]], base_name="$(nw)_sc_$(i)",
        start = PMD.comp_start_value(PMD.ref(pm, nw, :storage, i), ["sc_start", "sc"], c, 0.0)
    ) for i in PMD.ids(pm, nw, :storage))

    for (i, strg) in PMD.ref(pm, nw, :storage)
        if 4 in connections[i]
            ncnds = length(connections[i])-1
        else
            ncnds = length(connections[i])
        end 

        if bounded
            if haskey(strg, "charge_rating")
                for (idx, c) in enumerate(connections[i])
                    PMD.set_upper_bound(sc[i][c], max(strg["charge_rating"]/ncnds, 0.0))
                end
            end
            if haskey(strg, "discharge_rating")
                for (idx, c) in enumerate(connections[i])
                    PMD.set_lower_bound(sc[i][c], min(-1*strg["discharge_rating"]/ncnds, 0.0))
                end
            end
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :storage, :sc, PMD.ids(pm, nw, :storage), sc)
end

function optOn_variable_mc_generator_current(pm::PMD.AbstractUnbalancedIVRModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    optOn_variable_mc_generator_current_real(pm; nw=nw, bounded=bounded, report=report)
    optOn_variable_mc_generator_current_imaginary(pm; nw=nw, bounded=bounded, report=report)

    PMD.var(pm, nw)[:crg_bus] = Dict{Int, Any}()
    PMD.var(pm, nw)[:cig_bus] = Dict{Int, Any}()

    # store active and reactive power expressions for use in objective + post processing
    PMD.var(pm, nw)[:pg] = Dict{Int, Any}()
    PMD.var(pm, nw)[:qg] = Dict{Int, Any}()
end

"variable: `crg[j]` for `j` in `gen`"
function optOn_variable_mc_generator_current_real(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i,gen) in PMD.ref(pm, nw, :gen))
    crg = PMD.var(pm, nw)[:crg] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_crg_$(i)",
            start = PMD.comp_start_value(ref(pm, nw, :gen, i), "crg_start", c, 0.0)
        ) for i in PMD.ids(pm, nw, :gen)
    )
    if bounded
        for (i, g) in PMD.ref(pm, nw, :gen)
            if haskey(g, "crg")
                for (idx,c) in enumerate(connections[i])
                    JuMP.@constraint(pm.model, crg[i][c] == g["crg"][idx])
                end 
            else
                cmax = PMD._calc_gen_current_max(g, PMD.ref(pm, nw, :bus, g["gen_bus"]))
                for (idx,c) in enumerate(connections[i])
                    set_lower_bound(crg[i][c], -cmax[idx])
                    set_upper_bound(crg[i][c],  cmax[idx])
                end
            end 
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :gen, :crg, PMD.ids(pm, nw, :gen), crg)
end


"variable: `cig[j]` for `j` in `gen`"
function optOn_variable_mc_generator_current_imaginary(pm::PMD.AbstractUnbalancedPowerModel; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => gen["connections"] for (i,gen) in PMD.ref(pm, nw, :gen))
    cig = PMD.var(pm, nw)[:cig] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_cig_$(i)",
            start = PMD.comp_start_value(ref(pm, nw, :gen, i), "cig_start", c, 0.0)
        ) for i in PMD.ids(pm, nw, :gen)
    )
    if bounded
        for (i, g) in PMD.ref(pm, nw, :gen)
            if haskey(g, "cig")
                for (idx,c) in enumerate(connections[i])
                    JuMP.@constraint(pm.model, cig[i][c] == g["cig"][idx])
                end 
            else
                cmax = PMD._calc_gen_current_max(g, PMD.ref(pm, nw, :bus, g["gen_bus"]))
                for (idx,c) in enumerate(connections[i])
                    set_lower_bound(cig[i][c], -cmax[idx])
                    set_upper_bound(cig[i][c],  cmax[idx])
                end
            end 
        end
    end

    report && IM.sol_component_value(pm, PMD.pmd_it_sym, nw, :gen, :cig, PMD.ids(pm, nw, :gen), cig)
end