"""
    Kirchoff current law for buses.
    I_gen_r[j]: Real current injection from all generators connected at bus j.
                It is zero if no generator is connected at bus j.
    I_gen_i[j]: Real current injection from all generators connected at bus j.
                It is zero if no generator is connected at bus j.

"""


function kcl(Ybus, V_r, V_i, I_injections_r, I_injections_i)
    I_bus = Ybus*(V_r + V_i.*1im)
    return [I_injections_r - real(I_bus); I_injections_i - imag(I_bus)]
end
