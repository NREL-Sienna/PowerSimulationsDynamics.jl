
# TODO: Implement sparse methods when moving to large scale
struct SmallSignalOutput
    reduced_jacobian::Matrix{Float64}
    eigenvalues::Vector{Complex{Float64}}
    eigenvectors::Matrix{Complex{Float64}}
    index::Dict
    stable::Bool
    operating_point::Vector{Float64}
    damping::Dict{String, Dict{Symbol, Float64}}
    participation_factors::Dict{String, Dict{Symbol, Array{Float64}}}
end

function _determine_stability(vals::Vector{Complex{Float64}})
    for real_eig in real(vals)
        real_eig > 0.0 && return false
    end
    return true
end

function _calculate_forwardiff_jacobian(
    sim::Simulation{ImplicitModel},
    x_eval::Vector{Float64},
)
    var_count = get_variable_count(sim.simulation_inputs)
    dx0 = zeros(var_count) #Define a vector of zeros for the derivative
    sysf! = (out, x) -> system_implicit!(
        out,            #output of the function
        dx0,            #derivatives equal to zero
        x,              #states
        sim.simulation_inputs,     #Parameters
        0.0,            #time equals to zero.
    )
    out = zeros(var_count) #Define a vector of zeros for the output
    jacobian = ForwardDiff.jacobian(sysf!, out, x_eval)
    return jacobian
end

function _calculate_forwardiff_jacobian(
    sim::Simulation{MassMatrixModel},
    x_eval::Vector{Float64},
)
    var_count = get_variable_count(sim.simulation_inputs)
    sysf! = (dx, x) -> system_mass_matrix!(
        dx,            #derivatives equal to zero
        x,              #states
        sim.simulation_inputs,     #Parameters
        0.0,            #time equals to zero.
    )
    dx = zeros(var_count) #Define a vector of zeros for the output
    jacobian = ForwardDiff.jacobian(sysf!, dx, x_eval)
    return jacobian
end

function _make_reduce_jacobian_index(global_index, diff_states)
    jac_index = Dict{String, Dict{Symbol, Int}}()
    for (device_name, device_index) in global_index
        jac_index[device_name] = Dict{Symbol, Int}()
        for (state, ix) in device_index
            state_is_differential = diff_states[ix]
            if state_is_differential
                jac_index[device_name][state] = sum(diff_states[1:ix])
            elseif !state_is_differential
                #jac_index[device_name][state] = nothing
            else
                @assert false
            end
        end
    end
    return jac_index
end

function _get_state_types(sim::Simulation)
    #var_count = get_variable_count(sim.simulation_inputs)
    #bus_count = get_bus_count(sim.simulation_inputs)
    #diff_states = collect(trues(var_count))
    #diff_states[1:(2 * bus_count)] .= false
    #for b_ix in get_voltage_buses_ix(sim.simulation_inputs)
    #    diff_states[b_ix] = true
    #    diff_states[b_ix + bus_count] = true
    #end
    diff_states = copy(sim.simulation_inputs.DAE_vector)
    alg_states = .!diff_states
    jac_index =
        _make_reduce_jacobian_index(get_global_index(sim.simulation_inputs), diff_states)
    return jac_index, diff_states, alg_states
end

function _reduce_jacobian(
    jacobian::Matrix{Float64},
    diff_states::Vector{Bool},
    alg_states::BitArray{1},
    mass_matrix::SparseArrays.SparseMatrixCSC{Float64, Int},
)
    fx = @view jacobian[diff_states, diff_states]
    gy = jacobian[alg_states, alg_states]
    fy = @view jacobian[diff_states, alg_states]
    gx = @view jacobian[alg_states, diff_states]
    M = @view mass_matrix[diff_states, diff_states]
    inv_diag_M = Vector(1.0 ./ LinearAlgebra.diag(M))
    # TODO: Make operation using BLAS!
    reduced_jacobian = inv_diag_M .* (fx - fy * inv(gy) * gx)
    return reduced_jacobian
end

function _get_eigenvalues(reduced_jacobian::Matrix{Float64}, multimachine::Bool)
    eigen_vals, R_eigen_vect = LinearAlgebra.eigen(reduced_jacobian)
    if multimachine
        @warn(
            "No Infinite Bus found. Confirm stability directly checking eigenvalues.\nIf all eigenvalues are on the left-half plane and only one eigenvalue is zero, the system is small signal stable."
        )
        info_evals = "Eigenvalues are:\n"
        for i in eigen_vals
            info_evals = info_evals * string(i) * "\n"
        end
        @info(info_evals)
    end
    return eigen_vals, R_eigen_vect
end

function _get_damping(
    eigen_vals::Vector{Complex{Float64}},
    jac_index::Dict{String, Dict{Symbol, Int}},
)
    damping_results = Dict{String, Dict{Symbol, Float64}}()
    for (device_name, device_index) in jac_index
        damping_results[device_name] = Dict{Symbol, Float64}()
        for (state, ix) in device_index
            isnothing(ix) && continue
            eigen_val = eigen_vals[ix]
            damping_results[device_name][state] =
                -1 * real(eigen_val) / sqrt(real(eigen_val)^2 + imag(eigen_val)^2)
        end
    end
    return damping_results
end

function _get_participation_factors(
    R_eigen_vect::Matrix{Complex{Float64}},
    jac_index::Dict{String, Dict{Symbol, Int}},
)
    L_eigen_vect = inv(R_eigen_vect)
    participation_factors = Dict{String, Dict{Symbol, Array{Float64}}}()
    for (device_name, device_index) in jac_index
        participation_factors[device_name] = Dict{Symbol, Array{Float64}}()
        for (state, ix) in device_index
            den = sum(abs.(L_eigen_vect[:, ix]) .* abs.(R_eigen_vect[ix, :]))
            participation_factors[device_name][state] =
                abs.(L_eigen_vect[:, ix]) .* abs.(R_eigen_vect[ix, :]) ./ den
        end
    end
    return participation_factors
end

function small_signal_analysis(sim::Simulation; kwargs...)
    simulation_pre_step!(sim, get(kwargs, :reset_simulation, false), Real)
    sim.status = CONVERTED_FOR_SMALL_SIGNAL
    mass_matrix = get_mass_matrix(sim.simulation_inputs)
    x_eval = get(kwargs, :operating_point, sim.x0_init)
    jacobian = _calculate_forwardiff_jacobian(sim, x_eval)
    jac_index, diff_states, alg_states = _get_state_types(sim)
    reduced_jacobian = _reduce_jacobian(jacobian, diff_states, alg_states, mass_matrix)
    eigen_vals, R_eigen_vect = _get_eigenvalues(reduced_jacobian, sim.multimachine)
    damping = _get_damping(eigen_vals, jac_index)
    stable = _determine_stability(eigen_vals)
    participation_factors = _get_participation_factors(R_eigen_vect, jac_index)
    return SmallSignalOutput(
        reduced_jacobian,
        eigen_vals,
        R_eigen_vect,
        jac_index,
        stable,
        x_eval,
        damping,
        participation_factors,
    )
end
