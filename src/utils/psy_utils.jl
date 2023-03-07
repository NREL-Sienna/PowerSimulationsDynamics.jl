# This is where we do some type piracy on the PSY types
get_n_buses(sys::PSY.System) = length(sys.bus_numbers)

function _filter_function(x::T) where {T <: PSY.StaticInjection}
    if PSY.get_dynamic_injector(x) === nothing
        return false
    end

    if hasfield(T, :status)
        return PSY.get_status(x) * PSY.get_available(x)
    else
        return PSY.get_available(x)
    end
    error("Filtering of $(PSY.get_name(x)) failed")
end

function get_injectors_with_dynamics(sys::PSY.System)
    return PSY.get_components(x -> _filter_function(x), PSY.StaticInjection, sys)
end

function get_injection_without_dynamics(sys::PSY.System)
    return PSY.get_components(
        x ->
            PSY.get_dynamic_injector(x) === nothing &&
                PSY.get_available(x) &&
                !isa(x, PSY.ElectricLoad),
        PSY.StaticInjection,
        sys,
    )
end

function get_dynamic_branches(sys::PSY.System)
    return PSY.get_components(x -> PSY.get_available(x), PSY.DynamicBranch, sys)
end

function _transform_all_lines!(sys::PSY.System)
    for br in PSY.get_components(PSY.DynamicBranch, sys)
        dyn_br = DynamicBranch(br)
        @debug "Converted $(PSY.get_name(dyn_br)) to DynamicBranch"
        add_component!(sys, dyn_br)
    end
end

function transform_ybus_to_rectangular(
    ybus::SparseArrays.SparseMatrixCSC{Complex{Float64}, Int},
)
    # TODO: Improve performance here
    return hcat(vcat(real(ybus), -imag(ybus)), vcat(imag(ybus), real(ybus)))
end

function transform_branches_to_dynamic(sys::PSY.System, ::Type{T}) where {T <: PSY.ACBranch}
    for b in PSY.get_components(T, sys)
        dyn_branch = PSY.DynamicBranch(b)
        PSY.add_component!(sys, dyn_branch)
    end
    return
end
