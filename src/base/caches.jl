# The implementation of caches is influenced by the SparseDiffTools.jl. We implement this
# custom code since the structure of the system models is not compatible with the functionalities
# in SparseDiffTools

struct JacobianCache{T, U}
    du::Vector{T}
    du_dual::Vector{U}
    current_injections_r::Vector{T}
    current_injections_r_dual::Vector{U}
    current_injections_i::Vector{T}
    current_injections_current_dual::Vector{U}
    injection_ode::Vector{T}
    injection_ode_dual::Vector{U}
    branches_ode::Vector{T}
    branches_ode_dual::Vector{U}
    current_bus::Vector{T}
    current_bus_dual::Vector{U}
    current_balance::Vector{T}
    current_balance_dual::Vector{U}
    inner_vars::Vector{T}
    inner_vars::Vector{U}
end

get_du(jc::JacobianCache, ::Type{Float64}) = jc.du
get_current_injections_r(jc::JacobianCache, ::Type{Float64}) = jc.current_injections_r
get_current_injections_i(jc::JacobianCache, ::Type{Float64}) = jc._injections_i
get_injection_ode(jc::JacobianCache, ::Type{Float64}) = jc.injection_ode
get_branches_ode(jc::JacobianCache, ::Type{Float64}) = jc.banches_ode
get_current_bus(jc::JacobianCache, ::Type{Float64}) = jc.current_bus
get_current_balance(jc::JacobianCache, ::Type{Float64}) = jc.current_balance
get_inner_vars(jc::JacobianCache, ::Type{Float64}) = jc.inner_vars

get_du(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.du_dual)
get_current_injections_r(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.current_injections_r_dual)
get_current_injections_i(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.current_injections_i_dual)
get_injection_ode(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.cinjection_ode_dual)
get_branches_ode(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.cbanches_ode_dual)
get_current_bus(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.ccurrent_bus_dual)
get_current_balance(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.ccurrent_balance_dual)
get_inner_vars(jc::JacobianCache, ::Type{T}) where {T <: ForwardDiff.Dual} =
    reinterpret(T, jc.cinner_vars_dual)

struct SimCache{T}
    current_injections_r::Vector{T}
    current_injections_i::Vector{T}
    injection_ode::Vector{T}
    branches_ode::Vector{T}
    current_bus::Vector{T}
    current_balance::Vector{T}
    inner_vars::Vector{T}
end

get_du(sc::SimCache, ::Type{Float64}) = sc.du
get_current_injections_r(sc::SimCache, ::Type{Float64}) = sc.current_injections_r
get_current_injections_i(sc::SimCache, ::Type{Float64}) = sc._injections_i
get_injection_ode(sc::SimCache, ::Type{Float64}) = sc.injection_ode
get_branches_ode(sc::SimCache, ::Type{Float64}) = sc.banches_ode
get_current_bus(sc::SimCache, ::Type{Float64}) = sc.current_bus
get_current_balance(sc::SimCache, ::Type{Float64}) = sc.current_balance
get_inner_vars(sc::SimCache, ::Type{Float64}) = sc.inner_vars

function add_aux_arrays!(inputs::SimulationInputs, ::Type{T}) where {T <: Real}
    @debug "Auxiliary Arrays created with Type $(T)"
    bus_count = get_bus_count(inputs)
    get_aux_arrays(inputs)[1] = collect(zeros(T, bus_count))
    get_aux_arrays(inputs)[2] = collect(zeros(T, bus_count))
    get_aux_arrays(inputs)[3] = collect(zeros(T, get_n_injection_states(inputs)))
    get_aux_arrays(inputs)[4] = collect(zeros(T, get_n_branches_states(inputs)))
    get_aux_arrays(inputs)[5] = collect(zeros(Complex{T}, bus_count))
    get_aux_arrays(inputs)[6] = collect(zeros(T, 2 * bus_count))
    return
end
