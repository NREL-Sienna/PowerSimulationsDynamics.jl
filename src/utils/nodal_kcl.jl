out[algebraic_range] = kcl(sys.Ybus, V_r, V_i,
                           I_injections_r, I_injections_i)

for (ix, b) in enumerate(OMIB.buses)
   alg_curr = alg_currents(Y, V_r, V_i, ix)
   out_new[ix] = I_injections_r[ix] - alg_curr[1]
   out_new[ix+bus_size] = I_injections_i[ix] - alg_curr[2]
end

function alg_currents(Ybus, V_r, V_i, ix)
    current = transpose(V_r+V_i.*1im)*Ybus[ix,:]
    return real(current), imag(current)
end
