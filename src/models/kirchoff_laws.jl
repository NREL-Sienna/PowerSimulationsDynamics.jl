function Ybus_current_kirchoff(inputs, V_r, V_i, I_injections_r, I_injections_i)
    I_bus = get_aux_arrays(inputs)[5]
    I_balance = get_aux_arrays(inputs)[6]
    Ybus = get_Ybus(inputs) + get_total_shunts(inputs)
    # Note: BLAS doesn't work because the the type of Vr and Vi is not Matrix of Complex
    LinearAlgebra.mul!(I_bus, Ybus, (V_r + V_i .* 1im))
    bus_count = length(I_bus)
    for n in eachindex(I_bus)
        I_balance[n] = I_injections_r[n] - real(I_bus[n])
        I_balance[n + bus_count] = I_injections_i[n] - imag(I_bus[n])
    end

    return I_balance
end
