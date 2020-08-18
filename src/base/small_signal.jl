struct SmallSignalOutput
    reduced_jacobian::Matrix{Float64}
    eigenvalues::Vector{Complex{Float64}}
    eigenvectors::Matrix{Complex{Float64}}
    stable::Bool
    operating_point::Vector{Float64}
    damping::Dict{String, Dict{Symbol, Float64}}
end

function _calculate_forwardiff_jacobian(sim::Simulation, x_eval::Vector{Float64})
    system = get_system(sim.simulation_inputs)
    var_count = get_variable_count(sim.simulation_inputs)
    dx0 = zeros(var_count) #Define a vector of zeros for the derivative
    bus_count = get_bus_count(sim.simulation_inputs)
    sysf! = (out, x) -> system!(
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

function _reduce_jacobian(jacobian::Matrix{Float64}, sim::Simulation)
    var_count = get_variable_count(sim.simulation_inputs)
    diff_states = collect(trues(var_count))
    diff_states[1:(2 * bus_count)] .= false
    for b_ix in get_voltage_buses_ix(sim.simulation_inputs)
        diff_states[b_ix] = true
        diff_states[b_ix + bus_count] = true
    end
    alg_states = .!diff_states
    fx = @view jacobian[diff_states, diff_states]
    gy = jacobian[alg_states, alg_states]
    fy = @view jacobian[diff_states, alg_states]
    gx = @view jacobian[alg_states, diff_states]
    # TODO: Make operation using BLAS!
    reduced_jacobian = fx - fy * inv(gy) * gx
end

function _get_eigenvalues(reduced_jacobian::Matrix{Float64}, multimachine::bool)
    eigen_vals, R_eigen_vect = LinearAlgebra.eigen(reduced_jacobian)
    if multimachine
        @warn("No Infinite Bus found. Confirm stability directly checking eigenvalues.\nIf all eigenvalues are on the left-half plane and only one eigenvalue is zero, the system is small signal stable.")
        info_evals = "Eigenvalues are:\n"
        for i in eigen_vals
            info_evals = info_evals * string(i) * "\n"
        end
        @info(info_evals)
    end
    return eigen_vals, R_eigen_vect
end

function _get_damping(sim::Simulation, eigen_vals::Vector{Float64})
    damping_results = Dict{String, Dict{Symbol, Float64}}()
    for (device_name, device_index) in get_global_index(sim.simulation_inputs)
        damping_results[device_name] = Dict{Symbol, Float}()
        for (state, ix) in device_index
            eigen_val = eigen_vals[ix]
            damping_results[device_name][state] =
                -1 * real(eigen_val) / sqrt(real(eigen_val)^2 + imag(eigen_val)^2)
        end
    end
    return damping_results
end

#=
function _get_participation_factors(sim::Simulation, R_eigen_vect)
    participation_factors_results =
    L_eigen_vect = inv(R_eigen_vect)

    participation_factors = zeros(size(L_eigen_vect))
     den = sum(abs.(L_eigen_vect[:, ix]) .* abs.(R_eigen_vect[ix, :]))
    participation_factors[ix, :] =
            abs.(L_eigen_vect[:, ix]) .* abs.(R_eigen_vect[ix, :]) ./ den
        current_eigenvalue = real(eigen_val)

end
=#

function small_signal_analysis(sim::Simulation; kwargs...)
    reset_simulation = get(kwargs, :reset_simulation, false)
    _simulation_pre_step(sim, reset_simulation)
    _change_vector_type!(sim.simulation_inputs, Real)
    x_eval = get(kwargs, :operating_point, sim.x0_init)
    jacobian = _calculate_forwardiff_jacobian(sim, x_val)
    reduced_jacobian = _reduce_jacobian(jacobian, sim)
    eigen_vals, R_eigen_vect = _get_eigenvalues(reduced_jacobian, sim.multimachine)
    damping = _get_damping(sim, eigen_vals)
    stable = _determine_stability(eigen_vals)
    return SmallSignalOutput(
        reduced_jacobian,
        eigen_vals,
        R_eigen_vect,
        stable,
        x_eval,
        damping,
    )
end
