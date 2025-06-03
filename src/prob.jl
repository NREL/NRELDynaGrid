
function constraint_mc_storage_power_setpoint_imaginary(pm, i::Int)
    qs = [PMD.var(pm, :qs, i)[c] for c in PMD.ref(pm, :storage, i, "connections")]
    qs_ref = PMD.ref(pm, :storage, i, "qs")
    JuMP.@constraint(pm.model, sum(qs) == qs_ref)
end

function constraint_mc_gen_power_setpoint_imaginary(pm, i::Int)
    qg = [PMD.var(pm, :qg, i)[c] for c in PMD.ref(pm, :gen, i, "connections")]
    qg_ref = PMD.ref(pm, :gen, i, "qg")
    JuMP.@constraint(pm.model, qg .== qg_ref)
end


function build_mc_opf_oltc_dg_nc(pm::PMD.AbstractUnbalancedPowerModel)
    PMD.variable_mc_bus_voltage_on_off(pm; bounded=false)
    PMD.variable_mc_branch_power(pm; bounded=false)
    PMD.variable_mc_switch_power(pm; bounded=false)
    PMD.variable_mc_transformer_power(pm; bounded=false)
    PMD.variable_mc_generator_power_on_off(pm; bounded=false)
    PMD.variable_mc_load_power(pm; bounded=false)
    optOn_variable_mc_storage_power_mi_on_off(pm; bounded=false)
    
    # Need transformer taps
    PMD.variable_mc_oltc_transformer_tap(pm)
    PMD.variable_mc_capcontrol(pm; relax=true)
    
    PMD.constraint_mc_model_voltage(pm)
    
    for (i,bus) in PMD.ref(pm, :ref_buses)
    @assert bus["bus_type"] == 3
    
    PMD.constraint_mc_theta_ref(pm, i)
    PMD.constraint_mc_voltage_magnitude_only(pm, i)
    end
    
    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :gen)
    PMD.constraint_mc_generator_power(pm, id)
    end
    
    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :load)
    PMD.constraint_mc_load_power(pm, id)
    end
    
    for (i,bus) in PMD.ref(pm, :bus)
    PMD.constraint_mc_power_balance_capc(pm, i)
    end
    
    for i in PMD.ids(pm, :storage)
    PMD.constraint_storage_state(pm, i)
    PMD.constraint_storage_complementarity_nl(pm, i)
    PMD.constraint_mc_storage_losses(pm, i)
    PMD.constraint_mc_storage_thermal_limit(pm, i)
    end
    
    for i in PMD.ids(pm, :branch)
    PMD.constraint_mc_ohms_yt_from(pm, i)
    PMD.constraint_mc_ohms_yt_to(pm, i)
    end
    
    for i in PMD.ids(pm, :switch)
    PMD.constraint_mc_switch_state(pm, i)
    end
    
    for i in PMD.ids(pm, :transformer)
    PMD.constraint_mc_transformer_power(pm, i; fix_taps=false) # need to un-fix taps to get right voltages
    end
    
end
    
 
function build_mc_pf_oltc_dg_nc(pm::PMD.AbstractUnbalancedPowerModel)
    PMD.variable_mc_bus_voltage(pm; bounded=false)
    PMD.variable_mc_branch_power(pm; bounded=false)
    PMD.variable_mc_switch_power(pm; bounded=false)
    PMD.variable_mc_transformer_power(pm; bounded=false)
    PMD.variable_mc_generator_power(pm; bounded=false)
    PMD.variable_mc_load_power(pm; bounded=false)
    PMD.variable_mc_storage_power(pm; bounded=false)
    
    # Need transformer taps
    PMD.variable_mc_oltc_transformer_tap(pm)
    PMD.variable_mc_capcontrol(pm; relax=true)
    
    PMD.constraint_mc_model_voltage(pm)
    
    for (i,bus) in PMD.ref(pm, :ref_buses)
    @assert bus["bus_type"] == 3
    
    PMD.constraint_mc_theta_ref(pm, i)
    PMD.constraint_mc_voltage_magnitude_only(pm, i)
    end
    
    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :gen)
    PMD.constraint_mc_generator_power(pm, id)
    end
    
    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :load)
    PMD.constraint_mc_load_power(pm, id)
    end
    
    for (i,bus) in PMD.ref(pm, :bus)
    PMD.constraint_mc_power_balance_capc(pm, i)
    end
    
    for i in PMD.ids(pm, :storage)
    PMD.constraint_storage_state(pm, i)
    PMD.constraint_storage_complementarity_nl(pm, i)
    PMD.constraint_mc_storage_losses(pm, i)
    PMD.constraint_mc_storage_thermal_limit(pm, i)
    end
    
    for i in PMD.ids(pm, :branch)
    PMD.constraint_mc_ohms_yt_from(pm, i)
    PMD.constraint_mc_ohms_yt_to(pm, i)
    end
    
    for i in PMD.ids(pm, :switch)
    PMD.constraint_mc_switch_state(pm, i)
    end
    
    for i in PMD.ids(pm, :transformer)
    PMD.constraint_mc_transformer_power(pm, i; fix_taps=false) # need to un-fix taps to get right voltages
    end
end


function build_mc_pf_oltc_dg(pm::PMD.AbstractUnbalancedPowerModel)
    PMD.variable_mc_bus_voltage(pm; bounded=false)
    PMD.variable_mc_branch_power(pm; bounded=false)
    PMD.variable_mc_switch_power(pm; bounded=false)
    PMD.variable_mc_transformer_power(pm; bounded=false)
    PMD.variable_mc_generator_power(pm; bounded=true)
    PMD.variable_mc_load_power(pm; bounded=false)
    optOn_variable_mc_storage_power(pm; bounded=true)

    # Need transformer taps
    PMD.variable_mc_oltc_transformer_tap(pm)
    PMD.variable_mc_capcontrol(pm; relax=true)

    PMD.constraint_mc_model_voltage(pm)

    for (i,bus) in PMD.ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3

        PMD.constraint_mc_theta_ref(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
    end

    for (i,bus) in PMD.ref(pm, :bus)
        PMD.constraint_mc_power_balance_capc(pm, i)

        # PQ Bus Constraints
        if (length(PMD.ref(pm, :bus_gens, i)) > 0 || length(PMD.ref(pm, :bus_storages, i)) > 0) && !(i in PMD.ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            for j in PMD.ref(pm, :bus_gens, i)
                PMD.constraint_mc_gen_power_setpoint_real(pm, j)
                constraint_mc_gen_power_setpoint_imaginary(pm, j)
            end
            for j in PMD.ref(pm, :bus_storages, i)
                PMD.constraint_mc_storage_power_setpoint_real(pm, j)
                constraint_mc_storage_power_setpoint_imaginary(pm, j)
            end
        end
    end

    for i in PMD.ids(pm, :storage)
        optOn_constraint_storage_complementarity_nl(pm, i)
        optOn_constraint_mc_storage_losses(pm, i)
        PMD.constraint_mc_storage_thermal_limit(pm, i)
        optOn_constraint_storage_state(pm, i)
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_ohms_yt_from(pm, i)
        PMD.constraint_mc_ohms_yt_to(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_state(pm, i)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_power(pm, i; fix_taps=false)  # need to un-fix taps to get right voltages
    end

end

# Function identical to PMD.build_mc_opf_oltc except objective is changed to include storage cost
function build_mc_opf_oltc_dg(pm::PMD.AbstractUnbalancedPowerModel)
    PMD.variable_mc_bus_voltage(pm)

    PMD.variable_mc_branch_power(pm)
    PMD.variable_mc_switch_power(pm)
    PMD.variable_mc_transformer_power(pm)

    PMD.variable_mc_oltc_transformer_tap(pm)

    PMD.variable_mc_generator_power(pm)
    PMD.variable_mc_load_power(pm)
    PMD.variable_mc_storage_power(pm)

    PMD.constraint_mc_model_voltage(pm)

    for i in PMD.ids(pm, :ref_buses)
        PMD.constraint_mc_theta_ref(pm, i)
    end

    # generators should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_power_balance(pm, i)
    end

    for i in PMD.ids(pm, :storage)
        PMD.constraint_storage_state(pm, i)
        PMD.constraint_storage_complementarity_nl(pm, i)
        PMD.constraint_mc_storage_losses(pm, i)
        PMD.constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_ohms_yt_from(pm, i)
        PMD.constraint_mc_ohms_yt_to(pm, i)

        PMD.constraint_mc_voltage_angle_difference(pm, i)

        PMD.constraint_mc_thermal_limit_from(pm, i)
        PMD.constraint_mc_thermal_limit_to(pm, i)
        PMD.constraint_mc_ampacity_from(pm, i)
        PMD.constraint_mc_ampacity_to(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_state(pm, i)
        PMD.constraint_mc_switch_thermal_limit(pm, i)
        PMD.constraint_mc_switch_ampacity(pm, i)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_power(pm, i, fix_taps=false)
    end

    objective_mc_min_fuel_cost_w_storage(pm)
end

# # # # # # # # # # IVR form # # # # # # # # # # # # # #
"""
	function build_mc_opf(
		pm::AbstractUnbalancedIVRModel
	)

constructor for OPF in current-voltage variable space
"""
function build_mc_opf_oltc_dg(pm::PMD.AbstractUnbalancedIVRModel)
    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_switch_current(pm)
    optOn_variable_mc_transformer_current(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_load_current(pm)
    optOn_variable_mc_storage_power(pm)

    # Constraints
    for i in PMD.ids(pm, :ref_buses)
        PMD.constraint_mc_theta_ref(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)

        PMD.constraint_mc_bus_voltage_drop(pm, i)

        PMD.constraint_mc_voltage_angle_difference(pm, i)

        PMD.constraint_mc_thermal_limit_from(pm, i)
        PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_state(pm, i)
        PMD.constraint_mc_switch_current_limit(pm, i)
    end

    for i in PMD.ids(pm, :storage)
        optOn_constraint_mc_storage_losses(pm, i)
        optOn_constraint_storage_state(pm, i)
        optOn_constraint_storage_complementarity_nl(pm, i)
        PMD.constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_power(pm, i)
    end

    # Objective
    objective_mc_min_fuel_cost_w_storage(pm)
end

"""
	function build_mc_opf(
		pm::AbstractUnbalancedIVRModel
	)

constructor for OPF in current-voltage variable space
"""
function build_mc_pf_oltc_dg(pm::PMD.AbstractUnbalancedIVRModel)
    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_switch_current(pm)
    optOn_variable_mc_transformer_current(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_load_current(pm)
    optOn_variable_mc_storage_power(pm)

    # Constraints
    for i in PMD.ids(pm, :ref_buses)
        PMD.constraint_mc_theta_ref(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
    end


    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)

        PMD.constraint_mc_bus_voltage_drop(pm, i)

        PMD.constraint_mc_voltage_angle_difference(pm, i)

        PMD.constraint_mc_thermal_limit_from(pm, i)
        PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_state(pm, i)
        PMD.constraint_mc_switch_current_limit(pm, i)
    end

    for i in PMD.ids(pm, :storage)
        #PMD.constraint_mc_storage_thermal_limit(pm, i; nw=n)  #maybe add later 
        optOn_constraint_mc_storage_losses(pm, i)
        optOn_constraint_storage_state(pm, i)
        optOn_constraint_storage_complementarity_nl(pm, i)
        PMD.constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_power(pm, i)
    end

end

# # # # # # # # # # IVR form with ideal (lossless) storage # # # # # # # # # # # # # #
"""
	function build_mc_opf(
		pm::AbstractUnbalancedIVRModel
	)

constructor for OPF in current-voltage variable space
"""
function build_mc_opf_oltc_dg_idealStorage(pm::PMD.AbstractUnbalancedIVRModel)
    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_switch_current(pm)
    optOn_variable_mc_transformer_current(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_load_current(pm)
    optOn_variable_mc_storage_power_idealStorage(pm)


    # Constraints
    for i in PMD.ids(pm, :ref_buses)
        optOn_constraint_mc_theta_ref(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)

        PMD.constraint_mc_bus_voltage_drop(pm, i)

        PMD.constraint_mc_voltage_angle_difference(pm, i)

        PMD.constraint_mc_thermal_limit_from(pm, i)
        PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_state(pm, i)
        PMD.constraint_mc_switch_current_limit(pm, i)
    end

    for i in PMD.ids(pm, :storage)
        optOn_constraint_mc_storage_losses_idealStorage(pm, i)
        optOn_constraint_storage_state_idealStorage(pm, i)
        PMD.constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_power(pm, i)
    end

    # Objective
    objective_mc_min_fuel_cost_w_storage(pm)
end

"""
	function build_mc_opf(
		pm::AbstractUnbalancedIVRModel
	)

constructor for OPF in current-voltage variable space
"""
function build_mc_pf_oltc_dg_idealStorage(pm::PMD.AbstractUnbalancedIVRModel)
    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_switch_current(pm)
    optOn_variable_mc_transformer_current(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_load_current(pm)
    optOn_variable_mc_storage_power_idealStorage(pm)

    # Constraints
    for i in PMD.ids(pm, :ref_buses)
        optOn_constraint_mc_theta_ref(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)

        # PQ Bus Constraints
        if (length(PMD.ref(pm, :bus_gens, i)) > 0 || length(PMD.ref(pm, :bus_storages, i)) > 0) && !(i in PMD.ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            for j in PMD.ref(pm, :bus_gens, i)
                PMD.constraint_mc_gen_power_setpoint_real(pm, j)
                constraint_mc_gen_power_setpoint_imaginary(pm, j)
            end
            for j in PMD.ref(pm, :bus_storages, i)
                constraint_mc_storage_power_setpoint_real(pm, j)
                constraint_mc_storage_power_setpoint_imaginary(pm, j)
            end
            
        end
        #
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)

        PMD.constraint_mc_bus_voltage_drop(pm, i)

        PMD.constraint_mc_voltage_angle_difference(pm, i)

        PMD.constraint_mc_thermal_limit_from(pm, i)
        PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_state(pm, i)
        PMD.constraint_mc_switch_current_limit(pm, i)
    end

    for i in PMD.ids(pm, :storage)
        optOn_constraint_mc_storage_losses_idealStorage(pm, i)
        PMD.constraint_mc_storage_thermal_limit(pm, i)
        optOn_constraint_storage_state_idealStorage(pm, i)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_power(pm, i)
    end

    # Objective
    objective_mc_min_fuel_cost_w_storage(pm)
end
