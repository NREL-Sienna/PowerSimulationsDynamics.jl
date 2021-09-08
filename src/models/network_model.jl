function network_model(
    inputs::SimulationInputs,
    cache::Cache,
    voltages::AbstractArray{T},
) where {T <: Real}
    I_balance = get_current_balance(cache, T)
    # This operation might need improvement, when the total shunts aren't added the
    # function is not allocating
    ybus = get_ybus(inputs) .+ get_total_shunts(inputs)
    bus_count = get_bus_count(inputs)
    rows_vals = SparseArrays.rowvals(ybus)
    ybus_vals = SparseArrays.nonzeros(ybus)
    for n in UnitRange{Int}(1, bus_count)
        # The indexing works this way because of the shape of ybus_rectangular and because of the
        # ordering of
        for i in SparseArrays.nzrange(ybus, n + bus_count)
            I_balance[n] -= ybus_vals[i] * voltages[rows_vals[i]]
        end
        for i in SparseArrays.nzrange(ybus, n)
            I_balance[n + bus_count] -= ybus_vals[i] * voltages[rows_vals[i]]
        end
    end
    return I_balance
end

#=
A = sparse(I,J,V)
rows = rowvals(A)
vals = nonzeros(A)
m, n = size(A)
for j = 1:n
   for i in nzrange(A, j)
      row = rows[i]
      val = vals[i]
      # perform sparse wizardry...
   end
end
=#
