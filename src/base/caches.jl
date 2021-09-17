# The implementation of caches is influenced by the SparseDiffTools.jl. We implement this
# custom code since the structure of the system models is not compatible with the functionalities
# in SparseDiffTools

abstract type Cache end

struct JacobianCache{T, U} <: Cache
    bus_count::Int
    dx::Vector{T}
    dx_dual::Vector{U}
    injection_ode::Vector{T}
    injection_ode_dual::Vector{U}
    branches_ode::Vector{T}
    branches_ode_dual::Vector{U}
    current_balance::Vector{T}
    current_balance_dual::Vector{U}
    inner_vars::Vector{T}
    inner_vars_dual::Vector{U}
    global_vars::Vector{T}
    global_vars_dual::Vector{U}
end

function JacobianCache{T, U}(inputs::SimulationInputs) where {T, U}
    n = get_variable_count(inputs)
    n_inj = get_injection_n_states(inputs)
    n_branches = get_branches_n_states(inputs)
    bus_count = get_bus_count(inputs)
    inner_vars_count = get_inner_vars_count(inputs)
    n_global_vars = length(keys(get_global_vars_update_pointers(inputs)))
    return JacobianCache{T, U}(
        bus_count,
        zeros(T, n),
        zeros(U, n),
        zeros(T, n_inj),
        zeros(U, n_inj),
        zeros(T, n_branches),
        zeros(U, n_branches),
        zeros(T, 2 * bus_count),
        zeros(U, 2 * bus_count),
        zeros(T, inner_vars_count),
        zeros(U, inner_vars_count),
        zeros(T, n_global_vars),
        zeros(U, n_global_vars),
    )
end

get_dx(jc::JacobianCache, ::Type{Float64}) = jc.dx
get_current_injections_r(jc::JacobianCache, ::Type{Float64}) =
    jc.current_balance[1:(jc.bus_count)]
get_current_injections_i(jc::JacobianCache, ::Type{Float64}) =
    jc.current_balance[(jc.bus_count + 1):end]
get_injection_ode(jc::JacobianCache, ::Type{Float64}) = jc.injection_ode
get_branches_ode(jc::JacobianCache, ::Type{Float64}) = jc.branches_ode
get_current_balance(jc::JacobianCache, ::Type{Float64}) = jc.current_balance
get_inner_vars(jc::JacobianCache, ::Type{Float64}) = jc.inner_vars
get_global_vars(jc::JacobianCache, ::Type{Float64}) = jc.global_vars

get_du(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.dx_dual)
get_current_injections_r(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.current_balance_dual[1:(jc.bus_count)])
get_current_injections_i(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.current_balance_dual[(jc.bus_count + 1):end])
get_injection_ode(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.injection_ode_dual)
get_branches_ode(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.branches_ode_dual)
get_current_balance(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.current_balance_dual)
get_inner_vars(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.inner_vars_dual)
get_global_vars(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.global_vars_dual)

struct SimCache{T} <: Cache
    dx::Vector{T}
    current_injections_r::Vector{T}
    current_injections_i::Vector{T}
    injection_ode::Vector{T}
    branches_ode::Vector{T}
    current_balance::Vector{T}
    inner_vars::Vector{T}
    # always initialize with [1.0] for frequency reference value
    global_vars::Vector{T}
end

function SimCache{T}(inputs::SimulationInputs) where {T}
    n = get_variable_count(inputs)
    n_inj = get_injection_n_states(inputs)
    n_branches = get_branches_n_states(inputs)
    bus_count = get_bus_count(inputs)
    inner_vars_count = get_inner_vars_count(inputs)
    n_global_vars = length(keys(get_global_vars_update_pointers(inputs)))
    return SimCache{T}(
        zeros(T, n),
        zeros(T, bus_count),
        zeros(T, bus_count),
        zeros(T, n_inj),
        zeros(T, n_branches),
        zeros(T, bus_count),
        zeros(T, inner_vars_count),
        zeros(T, n_global_vars),
    )
end

get_dx(sc::SimCache, ::Type{Float64}) = sc.dx
get_current_injections_r(sc::SimCache, ::Type{Float64}) = sc.current_injections_r
get_current_injections_i(sc::SimCache, ::Type{Float64}) = sc._injections_i
get_injection_ode(sc::SimCache, ::Type{Float64}) = sc.injection_ode
get_branches_ode(sc::SimCache, ::Type{Float64}) = sc.branches_ode
get_current_balance(sc::SimCache, ::Type{Float64}) = sc.current_balance
get_inner_vars(sc::SimCache, ::Type{Float64}) = sc.inner_vars
get_global_vars(sc::SimCache, ::Type{Float64}) = sc.global_vars

get_ω_sys(cache::Cache, T::Type{<:Real}) =
    get_global_vars(cache, T)[GLOBAL_VAR_SYS_FREQ_INDEX]
