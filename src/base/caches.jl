# The implementation of caches is influenced by the SparseDiffTools.jl. We implement this
# custom code since the structure of the system models is not compatible with the functionalities
# in SparseDiffTools

abstract type Cache end

struct JacobianCache{F, U <: ForwardDiff.Dual} <: Cache
    f!::F
    bus_count::Int
    injection_ode::Vector{Float64}
    injection_ode_dual::Vector{U}
    branches_ode::Vector{Float64}
    branches_ode_dual::Vector{U}
    current_balance::Vector{Float64}
    current_balance_dual::Vector{U}
    inner_vars::Vector{Float64}
    inner_vars_dual::Vector{U}
    global_vars::Vector{Float64}
    global_vars_dual::Vector{U}
end

function JacobianCache{U}(F, inputs::SimulationInputs) where {U <: ForwardDiff.Dual}
    n = get_variable_count(inputs)
    n_inj = get_injection_n_states(inputs)
    n_branches = get_branches_n_states(inputs)
    bus_count = get_bus_count(inputs)
    inner_vars_count = get_inner_vars_count(inputs)
    n_global_vars = length(keys(get_global_vars_update_pointers(inputs)))
    return JacobianCache{typeof(F), U}(
        F,
        bus_count,
        zeros(Float64, n_inj),
        zeros(U, n_inj),
        zeros(Float64, n_branches),
        zeros(U, n_branches),
        zeros(Float64, 2 * bus_count),
        zeros(U, 2 * bus_count),
        zeros(Float64, inner_vars_count),
        zeros(U, inner_vars_count),
        setindex!(zeros(Float64, n_global_vars), 1.0, GLOBAL_VAR_SYS_FREQ_INDEX),
        setindex!(zeros(U, n_global_vars), convert(U, 1.0), GLOBAL_VAR_SYS_FREQ_INDEX),
    )
end

get_current_injections_r(jc::JacobianCache, ::Type{Float64}) =
    view(jc.current_balance, 1:(jc.bus_count))
get_current_injections_i(jc::JacobianCache, ::Type{Float64}) =
    view(jc.current_balance, (jc.bus_count + 1):(2 * jc.bus_count))
get_injection_ode(jc::JacobianCache, ::Type{Float64}) = jc.injection_ode
get_branches_ode(jc::JacobianCache, ::Type{Float64}) = jc.branches_ode
get_current_balance(jc::JacobianCache, ::Type{Float64}) = jc.current_balance
get_inner_vars(jc::JacobianCache, ::Type{Float64}) = jc.inner_vars
get_global_vars(jc::JacobianCache, ::Type{Float64}) = jc.global_vars

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

struct SimCache{F} <: Cache
    f!::F
    bus_count::Int
    injection_ode::Vector{Float64}
    branches_ode::Vector{Float64}
    current_balance::Vector{Float64}
    inner_vars::Vector{Float64}
    global_vars::Vector{Float64}
end

function SimCache(f!, inputs::SimulationInputs)
    n = get_variable_count(inputs)
    n_inj = get_injection_n_states(inputs)
    n_branches = get_branches_n_states(inputs)
    bus_count = get_bus_count(inputs)
    inner_vars_count = get_inner_vars_count(inputs)
    n_global_vars = length(keys(get_global_vars_update_pointers(inputs)))
    return SimCache{typeof(f!)}(
        f!,
        bus_count,
        zeros(Float64, n_inj),
        zeros(Float64, n_branches),
        zeros(Float64, 2 * bus_count),
        zeros(Float64, inner_vars_count),
        setindex!(zeros(Float64, n_global_vars), 1.0, GLOBAL_VAR_SYS_FREQ_INDEX),
    )
end

function get_current_injections_r(sc::SimCache, ::Type{Float64})
    return view(sc.current_balance, 1:(sc.bus_count))
end

function get_current_injections_i(sc::SimCache, ::Type{Float64})
    return view(sc.current_balance, ((sc.bus_count + 1):(sc.bus_count * 2)))
end

get_injection_ode(sc::SimCache, ::Type{Float64}) = sc.injection_ode
get_branches_ode(sc::SimCache, ::Type{Float64}) = sc.branches_ode
get_current_balance(sc::SimCache, ::Type{Float64}) = sc.current_balance
get_inner_vars(sc::SimCache, ::Type{Float64}) = sc.inner_vars
get_global_vars(sc::SimCache, ::Type{Float64}) = sc.global_vars

get_ω_sys(cache::Cache, T::Type{<:Union{Float64, ForwardDiff.Dual}}) =
    get_global_vars(cache, T)[GLOBAL_VAR_SYS_FREQ_INDEX]
