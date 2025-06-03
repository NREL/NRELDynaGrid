"""
    objective_mc_min_fuel_cost(pm::AbstractUnbalancedPowerModel)

Standard fuel cost minimization objective adapted to include cost for storage
"""
function objective_mc_min_fuel_cost_w_storage(pm::PMD.AbstractUnbalancedPowerModel; report::Bool=true)
    model = PMD.check_gen_cost_models(pm)

    if model == 1
        return objective_mc_min_fuel_cost_pwl_dg(pm; report=report)
    elseif model == 2
        return objective_mc_min_fuel_cost_polynomial_dg(pm; report=report)
    else
        error("Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end

"""
    objective_mc_min_fuel_cost_pwl(pm::AbstractUnbalancedPowerModel)

Fuel cost minimization objective with piecewise linear terms adapted to include cost for storage
"""
function objective_mc_min_fuel_cost_pwl_dg(pm::PMD.AbstractUnbalancedPowerModel; report::Bool=true)
    objective_mc_variable_pg_cost_dg(pm; report=report)
    return JuMP.@objective(pm.model, Min,
        sum(
            sum( PMD.var(pm, n, :pg_cost, i) for (i,gen) in nw_ref[:gen]) + sum( PMD.var(pm, n, :ps_cost, i) for (i,st) in nw_ref[:storage])
        for (n, nw_ref) in PMD.nws(pm))
    )
end

"""
    objective_mc_variable_pg_cost(pm::AbstractUnbalancedPowerModel)

adds pg_cost variables and constraints
"""
function objective_mc_variable_pg_cost_dg(pm::PMD.AbstractUnbalancedPowerModel; report::Bool=true)
    for (n, nw_ref) in PMD.nws(pm)
        pg_cost = PMD.var(pm, n)[:pg_cost] = Dict{Int,Any}()

        for (i,gen) in PMD.ref(pm, n, :gen)
            points = PMD.calc_pwl_points(gen["ncost"], gen["cost"], gen["pmin"], gen["pmax"])

            pg_cost_lambda = JuMP.@variable(pm.model,
                [i in 1:length(points)], base_name="$(n)_pg_cost_lambda",
                lower_bound = 0.0,
                upper_bound = 1.0
            )
            JuMP.@constraint(pm.model, sum(pg_cost_lambda) == 1.0)

            pg_expr = 0.0
            pg_cost_expr = 0.0
            for (i,point) in enumerate(points)
                pg_expr += point.mw*pg_cost_lambda[i]
                pg_cost_expr += point.cost*pg_cost_lambda[i]
            end
            JuMP.@constraint(pm.model, pg_expr == sum(PMD.var(pm, n, :pg, i)[c] for c in gen["connections"]))
            pg_cost[i] = pg_cost_expr
        end

        ps_cost = PMD.var(pm, n)[:ps_cost] = Dict{Int,Any}()

        for (i,st) in PMD.ref(pm, n, :storage)
            points = PMD.calc_pwl_points(st["ncost"], st["cost"], st["pmin"], st["pmax"])

            ps_cost_lambda = JuMP.@variable(pm.model,
                [i in 1:length(points)], base_name="$(n)_ps_cost_lambda",
                lower_bound = 0.0,
                upper_bound = 1.0
            )
            JuMP.@constraint(pm.model, sum(ps_cost_lambda) == 1.0)

            ps_expr = 0.0
            ps_cost_expr = 0.0
            for (i,point) in enumerate(points)
                ps_expr += point.mw*ps_cost_lambda[i]
                ps_cost_expr += point.cost*ps_cost_lambda[i]
            end
            JuMP.@constraint(pm.model, ps_expr == sum(PMD.var(pm, n, :ps, i)[c] for c in st["connections"]))
            ps_cost[i] = ps_cost_expr
        end

        report && _IM.sol_component_value(pm, pmd_it_sym, n, :gen, :pg_cost, PMD.ids(pm, n, :gen), pg_cost)
    end
end



"""
    objective_mc_min_fuel_cost_polynomial(pm::AbstractUnbalancedPowerModel)

Fuel cost minimization objective for polynomial terms adapted to include cost for storage
"""
function objective_mc_min_fuel_cost_polynomial_dg(pm::PMD.AbstractUnbalancedPowerModel; report::Bool=true)
    order = PMD.calc_max_cost_index(pm.data)-1

    if order <= 2
        return _objective_mc_min_fuel_cost_polynomial_linquad_dg(pm; report=report)
    else
        return _objective_mc_min_fuel_cost_polynomial_nl_dg(pm; report=report)
    end
end

"gen connections adaptation of min fuel cost polynomial linquad objective adapted to include cost for storage"
function _objective_mc_min_fuel_cost_polynomial_linquad_dg(pm::PMD.AbstractUnbalancedPowerModel; report::Bool=true)
    pg_contains_nl_exp = any(x<:JuMP.NonlinearExpression for x in vcat([typeof.(isa(pg, JuMP.Containers.DenseAxisArray) ? pg.data : pg) for nw in PMD.nw_ids(pm) for (id,pg) in PMD.var(pm, nw, :pg)]...))
    gen_cost = Dict()

    ps_contains_nl_exp = any(x<:JuMP.NonlinearExpression for x in vcat([typeof.(isa(ps, JuMP.Containers.DenseAxisArray) ? ps.data : ps) for nw in PMD.nw_ids(pm) for (id,ps) in PMD.var(pm, nw, :ps)]...))
    st_cost = Dict()

    if !pg_contains_nl_exp && !ps_contains_nl_exp
        for (n, nw_ref) in PMD.nws(pm)
            for (i,gen) in nw_ref[:gen]
                pg = sum(PMD.var(pm, n, :pg, i))

                if length(gen["cost"]) == 1
                    gen_cost[(n,i)] = gen["cost"][1]
                elseif length(gen["cost"]) == 2
                    gen_cost[(n,i)] = gen["cost"][1]*pg + gen["cost"][2]
                elseif length(gen["cost"]) == 3
                    gen_cost[(n,i)] = gen["cost"][1]*pg^2 + gen["cost"][2]*pg + gen["cost"][3]
                else
                    gen_cost[(n,i)] = 0.0
                end
            end
        end

        for (n, nw_ref) in PMD.nws(pm)
            for (i,st) in nw_ref[:storage]
                ps = sum(PMD.var(pm, n, :ps, i))

                if length(st["cost"]) == 1
                    st_cost[(n,i)] = st["cost"][1]
                elseif length(st["cost"]) == 2
                    st_cost[(n,i)] = st["cost"][1]*ps + st["cost"][2]
                elseif length(st["cost"]) == 3
                    st_cost[(n,i)] = st["cost"][1]*ps^2 + st["cost"][2]*ps + st["cost"][3]
                else
                    st_cost[(n,i)] = 0.0
                end
            end
        end

        return JuMP.@objective(pm.model, Min,
            sum(
                sum( gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] ) + sum( st_cost[(n,i)] for (i,st) in nw_ref[:storage] )
            for (n, nw_ref) in PMD.nws(pm))
        )
    else
        for (n, nw_ref) in PMD.nws(pm)
            for (i,gen) in nw_ref[:gen]
                bus = gen["gen_bus"]

                #to avoid function calls inside of @NLconstraint:
                pg = PMD.var(pm, n, :pg, i)
                pg = isa(pg, JuMP.Containers.DenseAxisArray) ? pg.data : pg

                int_dim = length(pg)
                if length(gen["cost"]) == 1
                    gen_cost[(n,i)] = gen["cost"][1]
                elseif length(gen["cost"]) == 2
                    gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*sum(pg[i] for i in 1:int_dim) + gen["cost"][2])
                elseif length(gen["cost"]) == 3
                    gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, gen["cost"][1]*sum(pg[i] for i in 1:int_dim)^2 + gen["cost"][2]*sum(pg[i] for i in 1:int_dim) + gen["cost"][3])
                else
                    gen_cost[(n,i)] = 0.0
                end
            end

            for (i,st) in nw_ref[:storage]
                bus = st["storage_bus"]

                #to avoid function calls inside of @NLconstraint:
                ps = PMD.var(pm, n, :ps, i)
                ps = isa(ps, JuMP.Containers.DenseAxisArray) ? ps.data : ps

                int_dim = length(ps)
                if length(st["cost"]) == 1
                    st_cost[(n,i)] = st["cost"][1]
                elseif length(st["cost"]) == 2
                    st_cost[(n,i)] = JuMP.@NLexpression(pm.model, st["cost"][1]*sum(ps[i] for i in 1:int_dim) + st["cost"][2])
                elseif length(st["cost"]) == 3
                    st_cost[(n,i)] = JuMP.@NLexpression(pm.model, st["cost"][1]*sum(1000*ps[i] for i in 1:int_dim)^2 + st["cost"][2]*sum(ps[i] for i in 1:int_dim) + st["cost"][3])
                else
                    st_cost[(n,i)] = 0.0
                end
            end
        end

        return JuMP.@NLobjective(pm.model, Min,
            sum(
                sum(    gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] ) + sum( st_cost[(n,i)] for (i,st) in nw_ref[:storage] )
            for (n, nw_ref) in PMD.nws(pm))
        )
    end
end

function _objective_mc_min_fuel_cost_polynomial_nl_dg(pm::PMD.AbstractUnbalancedPowerModel; report::Bool=true)
    gen_cost = Dict()
    st_cost = Dict()
    for (n, nw_ref) in PMD.nws(pm)
        for (i,gen) in nw_ref[:gen]
            pg = sum( PMD.var(pm, n, :pg, i)[c] for c in gen["connections"] )

            cost_rev = reverse(gen["cost"])
            if length(cost_rev) == 1
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1])
            elseif length(cost_rev) == 2
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*pg)
            elseif length(cost_rev) == 3
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*pg + cost_rev[3]*pg^2)
            elseif length(cost_rev) >= 4
                cost_rev_nl = cost_rev[4:end]
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*pg + cost_rev[3]*pg^2 + sum( v*pg^(d+2) for (d,v) in enumerate(cost_rev_nl)) )
            else
                gen_cost[(n,i)] = JuMP.@NLexpression(pm.model, 0.0)
            end
        end

        for (i,st) in nw_ref[:storage]
            ps = sum( PMD.var(pm, n, :ps, i)[c] for c in st["connections"] )

            cost_rev = reverse(st["cost"])
            if length(cost_rev) == 1
                st_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1])
            elseif length(cost_rev) == 2
                st_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*ps)
            elseif length(cost_rev) == 3
                st_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*ps + cost_rev[3]*ps^2)
            elseif length(cost_rev) >= 4
                cost_rev_nl = cost_rev[4:end]
                st_cost[(n,i)] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*ps + cost_rev[3]*ps^2 + sum( v*ps^(d+2) for (d,v) in enumerate(cost_rev_nl)) )
            else
                st_cost[(n,i)] = JuMP.@NLexpression(pm.model, 0.0)
            end
        end
    end

    return JuMP.@NLobjective(pm.model, Min,
        sum(
            sum( gen_cost[(n,i)] for (i,gen) in nw_ref[:gen] ) + sum( st_cost[(n,i)] for (i,st) in nw_ref[:storage] )
        for (n, nw_ref) in PMD.nws(pm))
    )
end

