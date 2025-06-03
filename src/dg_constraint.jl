

""
function optOn_constraint_storage_state_initial(pm::PMD.AbstractUnbalancedPowerModel, n::Int, i::Int, energy::Real, charge_eff::Real, discharge_eff::Real, time_elapsed::Real, connections::Vector{Int})
    sc = PMD.var(pm, n, :sc, i)
    sd = PMD.var(pm, n, :sd, i)
    se = PMD.var(pm, n, :se, i)

    JuMP.@constraint(pm.model, se - energy == time_elapsed*(charge_eff*sum(sc[c] for c in connections) - sum(sd[c] for c in connections)/discharge_eff))
    nothing
end

""
function optOn_constraint_storage_state(pm::PMD.AbstractUnbalancedPowerModel, n_1::Int, n_2::Int, i::Int, charge_eff::Real, discharge_eff::Real, time_elapsed::Real, connections::Vector{Int})
    sc_2 = var(pm, n_2, :sc, i)
    sd_2 = var(pm, n_2, :sd, i)
    se_2 = var(pm, n_2, :se, i)
    se_1 = var(pm, n_1, :se, i)

    JuMP.@constraint(pm.model, se_2 - se_1 == time_elapsed*(charge_eff*sum(sc_2[c] for c in connections) - sum(sd_2[c] for c in connections)/discharge_eff))
    nothing
end

"""
    constraint_storage_complementarity_nl(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)

Template function for nonlinear storage complementarity constraints
"""
function optOn_constraint_storage_complementarity_nl(pm::PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    storage = PMD.ref(pm, nw, :storage, i)
    optOn_constraint_storage_complementarity_nl(pm, nw, i, storage["connections"])
    nothing
end

""
function optOn_constraint_storage_complementarity_nl(pm::PMD.AbstractUnbalancedPowerModel, n::Int, i::Int, connections::Vector{Int})
    sc = var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)

    JuMP.@constraint(pm.model, sum(sc[c] for c in connections)*sum(sd[c] for c in connections) == 0.0)
    nothing
end

""
function optOn_constraint_mc_storage_losses(pm::PMD.AbstractUnbalancedPowerModel, nw::Int, i::Int, bus::Int, connections::Vector{Int}, r::Real, x::Real, p_loss::Real, q_loss::Real)
    storage = PMD.ref(pm, nw, :storage, i)

    vr  = PMD.var(pm, nw,  :vr, bus)
    vi  = PMD.var(pm, nw,  :vi, bus)
    ps  = PMD.var(pm, nw,  :ps, i)
    qs  = PMD.var(pm, nw,  :qs, i)
    sc  = PMD.var(pm, nw,  :sc, i)
    sd  = PMD.var(pm, nw,  :sd, i)
    qsc = PMD.var(pm, nw, :qsc, i)

    JuMP.@NLconstraint(pm.model,
        sum(ps[c] for c in connections) + sum((sd[c] - sc[c]) for c in connections)
        ==
        p_loss + r * sum((ps[c]^2 + qs[c]^2)/(vr[c]^2 + vi[c]^2) for (idx,c) in enumerate(connections))
    )

    JuMP.@NLconstraint(pm.model,
        sum(qs[c] for c in connections)
        ==
        sum(qsc[c] for c in connections) + q_loss + x * sum((ps[c]^2 + qs[c]^2)/(vr[c]^2 + vi[c]^2) for (idx,c) in enumerate(connections))
    )
end



# # # # # # # # # # Constraints for IVR form # # # # # # # # # # # # # #

"Creates phase angle constraints at reference buses"
function optOn_constraint_mc_theta_ref(pm::PMD.AbstractUnbalancedIVRModel, nw::Int, i::Int, va_ref::Vector{<:Real})
    terminals = PMD.ref(pm, nw, :bus, i, "terminals")
    vr = PMD.var(pm, nw, :vr, i)
    vi = PMD.var(pm, nw, :vi, i)

    # deal with cases first where tan(theta)==Inf or tan(theta)==0
    for (idx, t) in enumerate(terminals)
        if va_ref[idx] == pi/2
            JuMP.@constraint(pm.model, vr[t] == 0)
            JuMP.@constraint(pm.model, vi[t] >= 0)
        elseif va_ref[idx] == -pi/2
            JuMP.@constraint(pm.model, vr[t] == 0)
            JuMP.@constraint(pm.model, vi[t] <= 0)
        elseif va_ref[idx] == 0
            JuMP.@constraint(pm.model, vr[t] >= 0)
            JuMP.@constraint(pm.model, vi[t] == 0)
        elseif va_ref[idx] == pi
            JuMP.@constraint(pm.model, vr[t] >= 0)
            JuMP.@constraint(pm.model, vi[t] == 0)
        else
            JuMP.@constraint(pm.model, vi[t] == tan(va_ref[idx])*vr[t])
            # va_ref also implies a sign for vr, vi
            if 0<=va_ref[idx] && va_ref[idx] <= pi
                JuMP.@constraint(pm.model, vi[t] >= 0)
            else
                JuMP.@constraint(pm.model, vi[t] <= 0)
            end
        end
    end
end

function optOn_constraint_mc_storage_losses(pm::PMD.AbstractUnbalancedIVRModel, n::Int, i::Int, bus::Int, connections::Vector{Int}, r::Real, x::Real, p_loss::Real, q_loss::Real)
    
    if 4 in connections
        ncnds = length(connections)-1
    else
        ncnds = length(connections)
    end
    
    vr = PMD.var(pm, n, :vr, bus)
    vi = PMD.var(pm, n, :vi, bus)
    cr = PMD.var(pm, n, :crs, i)
    ci = PMD.var(pm, n, :cis, i)
    sc = PMD.var(pm, n, :sc, i)
    sd = var(pm, n, :sd, i)
    qsc = PMD.var(pm, n, :qsc, i)
    ps  = PMD.var(pm, n,  :ps, i)
    qs  = PMD.var(pm, n,  :qs, i)

    for c in connections
        JuMP.@constraint(pm.model, (vr[c]*cr[c] + vi[c]*ci[c]) == ps[c])
        JuMP.@constraint(pm.model, (vi[c]*cr[c] - vr[c]*ci[c]) == qs[c])

        JuMP.@constraint(pm.model, (vr[c]*cr[c] + vi[c]*ci[c]) + (sd[c] - sc[c])
        == (p_loss/ncnds) + r * ((cr[c]^2 + ci[c]^2)))

        JuMP.@constraint(pm.model, (vi[c]*cr[c] - vr[c]*ci[c]) 
        == qsc[c] + (q_loss/ncnds) + x * (cr[c]^2 + ci[c]^2)
    )
    end
end

@doc raw"""
    constraint_mc_storage_power_setpoint_real(pm::AbstractUnbalancedPowerModel, nw::Int, i::Int, ps::Real)::Nothing

Generic storage real power setpoint constraint

```math
P_s == P_s^{setpoint}
```
"""
function constraint_mc_storage_power_setpoint_real(pm::AbstractUnbalancedIVRModel, nw::Int, i::Int, ps::Vector{Float64})
    ps_var = [var(pm, nw, :ps, i)[c] for c in ref(pm, nw, :storage, i)["connections"]]

    JuMP.@constraint(pm.model, ps_var .== ps)
end

function constraint_mc_storage_power_setpoint_imaginary(pm::AbstractUnbalancedIVRModel, nw::Int, i::Int, qs::Vector{Float64})
    qs_var = [var(pm, nw, :qs, i)[c] for c in ref(pm, nw, :storage, i)["connections"]]

    JuMP.@constraint(pm.model, qs_var .== qs)
end

"`qq[i] == qq`"
function constraint_mc_gen_power_setpoint_imaginary(pm::AbstractUnbalancedIVRModel, n::Int, i, qg_ref)
    gen = ref(pm, n, :gen, i)
    bus = gen["gen_bus"]
    connections = gen["connections"]
    vr = [var(pm, n, :vr, bus)[c] for c in connections]
    vi = [var(pm, n, :vi, bus)[c] for c in connections]
    cr = [var(pm, n, :crg, i)[c] for c in connections]
    ci = [var(pm, n, :cig, i)[c] for c in connections]

    JuMP.@constraint(pm.model, qg_ref .== vi.*cr  - vr.*ci)
end 

function constraint_mc_gen_current_setpoint_real(pm::AbstractUnbalancedIVRModel, n::Int, i, crg_ref)
    gen = ref(pm, n, :gen, i)
    bus = gen["gen_bus"]
    connections = gen["connections"]
    cr = [var(pm, n, :crg, i)[c] for c in connections]

    JuMP.@constraint(pm.model, crg_ref .== cr )
end 

function constraint_mc_gen_current_setpoint_imaginary(pm::AbstractUnbalancedIVRModel, n::Int, i, cig_ref)
    gen = ref(pm, n, :gen, i)
    bus = gen["gen_bus"]
    connections = gen["connections"]
    ci = [var(pm, n, :cig, i)[c] for c in connections]

    JuMP.@constraint(pm.model, cig_ref .== ci)
end 

# # # # # # # # # # for IVR form with ideal (lossless) storage # # # # # # # # # # # # # #
function optOn_constraint_mc_storage_losses_idealStorage(pm::PMD.AbstractUnbalancedIVRModel, n::Int, i::Int, bus::Int, connections::Vector{Int}, r::Real, x::Real, p_loss::Real, q_loss::Real)
    
    if 4 in connections
        ncnds = length(connections)-1
    else
        ncnds = length(connections)
    end
    
    vr = PMD.var(pm, n, :vr, bus)
    vi = PMD.var(pm, n, :vi, bus)
    cr = PMD.var(pm, n, :crs, i)
    ci = PMD.var(pm, n, :cis, i)
    sc = PMD.var(pm, n, :sc, i)
    qsc = PMD.var(pm, n, :qsc, i)
    ps  = PMD.var(pm, n,  :ps, i)
    qs  = PMD.var(pm, n,  :qs, i)

    for c in connections
        JuMP.@constraint(pm.model, (vr[c]*cr[c] + vi[c]*ci[c]) == ps[c])
        JuMP.@constraint(pm.model, (vi[c]*cr[c] - vr[c]*ci[c]) == qs[c])

        JuMP.@constraint(pm.model, (vr[c]*cr[c] + vi[c]*ci[c]) - sc[c]
        == (p_loss/ncnds) + r * ((cr[c]^2 + ci[c]^2)))

        JuMP.@constraint(pm.model, (vi[c]*cr[c] - vr[c]*ci[c]) 
        == qsc[c] + (q_loss/ncnds) + x * (cr[c]^2 + ci[c]^2)
    )
    end

end


""
function optOn_constraint_storage_state_initial_idealStorage(pm::PMD.AbstractUnbalancedPowerModel, n::Int, i::Int, energy::Real, time_elapsed::Real, connections::Vector{Int})
    sc = PMD.var(pm, n, :sc, i)
    se = PMD.var(pm, n, :se, i)

    JuMP.@constraint(pm.model, se - energy == time_elapsed*sum(sc[c] for c in connections))
    nothing
end

""
function optOn_constraint_storage_state_idealStorage(pm::PMD.AbstractUnbalancedPowerModel, n_1::Int, n_2::Int, i::Int, time_elapsed::Real, connections::Vector{Int})
    sc_2 = var(pm, n_2, :sc, i)
    se_2 = var(pm, n_2, :se, i)
    se_1 = var(pm, n_1, :se, i)

    JuMP.@constraint(pm.model, se_2 - se_1 == time_elapsed*sum(sc_2[c] for c in connections))
    nothing
end