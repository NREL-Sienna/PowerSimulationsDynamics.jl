function limit_output_current(limiter :: Nothing, Id_cnv_ref :: Union{Float64, ForwardDiff.Dual}, Iq_cnv_ref :: Union{Float64, ForwardDiff.Dual})
    return Id_cnv_ref, Iq_cnv_ref
end

function limit_output_current(limiter:: PSY.InstantaneousOutputCurrentLimiter, Id_cnv_ref :: Union{Float64, ForwardDiff.Dual}, Iq_cnv_ref :: Union{Float64, ForwardDiff.Dual})
    d_lim = PSY.get_Id_max(limiter)
    q_lim = PSY.get_Iq_max(limiter)
    Id_cnv_ref2 = clamp(Id_cnv_ref, -d_lim, d_lim)
    Iq_cnv_ref2 = clamp(Iq_cnv_ref, -q_lim, q_lim)
    return Id_cnv_ref2, Iq_cnv_ref2
end

function limit_output_current(limiter :: PSY.MagnitudeOutputCurrentLimiter, Id_cnv_ref :: Union{Float64, ForwardDiff.Dual}, Iq_cnv_ref :: Union{Float64, ForwardDiff.Dual})
    limit_value = PSY.get_I_max(limiter)
    if (Id_cnv_ref^2 + Iq_cnv_ref^2)^(1/2) > limit_value
        Id_cnv_ref2 = limit_value*cos(theta)
        Iq_cnv_ref2 = limit_value*sin(theta)
    else
        Id_cnv_ref2 = Id_cnv_ref
        Iq_cnv_ref2 = Iq_cnv_ref
    end
    return Id_cnv_ref2, Iq_cnv_ref2
end

function limit_output_current(limiter :: PSY.SaturationOutputCurrentLimiter, Id_cnv_ref :: Union{Float64, ForwardDiff.Dual}, Iq_cnv_ref :: Union{Float64, ForwardDiff.Dual})
    limit_value = PSY.get_I_max(limiter)
    gain = PSY.get_kw(limiter)
    theta = atan(Iq_cnv_ref, Id_cnv_ref)
    if (Id_cnv_ref^2 + Iq_cnv_ref^2)^(1/2) > limit_value 
        Id_cnv_ref2 = limit_value*cos(theta)
        Iq_cnv_ref2 = limit_value*sin(theta)
    else
        Id_cnv_ref2 = Id_cnv_ref
        Iq_cnv_ref2 = Iq_cnv_ref
    end
    Del_Vv_d = gain * (Id_cnv_ref - Id_cnv_ref2)
    Del_Vv_q = gain * (Iq_cnv_ref - Iq_cnv_ref2)
    return Id_cnv_ref2, Iq_cnv_ref2, Del_Vv_d, Del_Vv_q
end

function limit_output_current(limiter :: PSY.HybridOutputCurrentLimiter, Id_cnv_ref :: Union{Float64, ForwardDiff.Dual}, Iq_cnv_ref :: Union{Float64, ForwardDiff.Dual}, ω :: Union{Float64, ForwardDiff.Dual})
    limit_value = PSY.get_I_max(limiter)
    real_imped = PSY.get_rv(limiter)
    imag_imped = PSY.get_lv(limiter)
    theta = atan(Iq_cnv_ref, Id_cnv_ref)
    if (Id_cnv_ref^2 + Iq_cnv_ref^2)^(1/2) > limit_value 
        Id_cnv_ref2 = limit_value*cos(theta)
        Iq_cnv_ref2 = limit_value*sin(theta)
    else
        Id_cnv_ref2 = Id_cnv_ref
        Iq_cnv_ref2 = Iq_cnv_ref
    end
    Del_Vv_d = ω * imag_imped * (Iq_cnv_ref2 - Iq_cnv_ref) + real_imped * (Id_cnv_ref - Id_cnv_ref2)
    Del_Vv_q = ω * imag_imped * (Id_cnv_ref - Id_cnv_ref2) + real_imped * (Iq_cnv_ref - Iq_cnv_ref2)
    return Id_cnv_ref2, Iq_cnv_ref2, Del_Vv_d, Del_Vv_q
end