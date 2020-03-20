"""
    Kirchoff law for buses.
    I_gen_r[j]: Real current injection from all generators connected at bus j.
                It is zero if no generator is connected at bus j.
    I_gen_i[j]: Real current injection from all generators connected at bus j.
                It is zero if no generator is connected at bus j.

"""
#TODO: Improve performance of this function
function kirchoff_laws(sys, V_r, V_i, I_injections_r, I_injections_i)
    Ybus = PSY.get_ext(sys)[YBUS]
    I_bus = Ybus * (V_r + V_i .* 1im)
    I_balance = [real(I_bus) - I_injections_r; imag(I_bus) - I_injections_i]
    return I_balance
end


#=
function kvl(I_net, B_eq)
    return [I_net[1].*B_eq ;
            I_net[2].*B_eq]
end
=#
