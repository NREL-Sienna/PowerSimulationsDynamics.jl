# The implementation of caches is influenced by the SparseDiffTools.jl. We implement this
# custom code since the structure of the system models is not compatible with the functionalities
# in SparseDiffTools

abstract type Cache end

struct JacobianCache{F, U <: ForwardDiff.Dual} <: Cache
    f!::F
    bus_count::Int
    ode_output::Vector{Float64}
    ode_output_dual::Vector{U}
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
    n_inj = get_injection_n_states(inputs)
    n_branches = get_branches_n_states(inputs)
    @debug "injection ode size = $n_inj branches ode size = $n_branches"
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
get_ode_output(jc::JacobianCache, ::Type{Float64}) = jc.ode_output
get_branches_ode(jc::JacobianCache, ::Type{Float64}) = jc.branches_ode
get_current_balance(jc::JacobianCache, ::Type{Float64}) = jc.current_balance
get_inner_vars(jc::JacobianCache, ::Type{Float64}) = jc.inner_vars
get_global_vars(jc::JacobianCache, ::Type{Float64}) = jc.global_vars

get_current_injections_r(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.current_balance_dual[1:(jc.bus_count)])
get_current_injections_i(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.current_balance_dual[(jc.bus_count + 1):end])
get_ode_output(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.ode_output_dual)
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
    ode_output::Vector{Real}
    branches_ode::Vector{Real}
    current_balance::Vector{Real}
    inner_vars::Vector{Real}
    global_vars::Vector{Real}
end

function SimCache(f!, inputs::SimulationInputs)
    n_inj = get_injection_n_states(inputs)
    n_branches = get_branches_n_states(inputs)
    @debug "injection ode size = $n_inj branches ode size = $n_branches"
    bus_count = get_bus_count(inputs)
    inner_vars_count = get_inner_vars_count(inputs)
    n_global_vars = length(keys(get_global_vars_update_pointers(inputs)))
    global_vars = setindex!(
        zeros(Real, n_global_vars),
        1.0,
        GLOBAL_VAR_SYS_FREQ_INDEX,
    )
    return SimCache{typeof(f!)}(
        f!,
        bus_count,
        zeros(Real, n_inj),
        zeros(Real, n_branches),
        zeros(Real, 2 * bus_count),
        zeros(Real, inner_vars_count),
        global_vars,
    )
end

function get_current_injections_r(sc::SimCache, ::Type{Float64})
    return view(sc.current_balance, 1:(sc.bus_count))
end

function get_current_injections_i(sc::SimCache, ::Type{Float64})
    return view(sc.current_balance, ((sc.bus_count + 1):(sc.bus_count * 2)))
end

get_ode_output(sc::SimCache, ::Type{T}) where {T <: Real} = sc.ode_output
get_branches_ode(sc::SimCache, ::Type{T}) where {T <: Real} = sc.branches_ode
get_current_balance(sc::SimCache, ::Type{T}) where {T <: Real} = sc.current_balance
get_inner_vars(sc::SimCache, ::Type{T}) where {T <: Real} = sc.inner_vars
get_global_vars(sc::SimCache, ::Type{T}) where {T <: Real} = sc.global_vars

get_ω_sys(cache::Cache, T::Type{<:ACCEPTED_REAL_TYPES}) =
    get_global_vars(cache, T)[GLOBAL_VAR_SYS_FREQ_INDEX]
