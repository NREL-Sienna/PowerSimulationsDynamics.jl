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
    theta = atan(Iq_cnv_ref, Id_cnv_ref)
    if (Id_cnv_ref^2 + Iq_cnv_ref^2)^(1/2) > limit_value 
        Id_cnv_ref2 = limit_value*cos(theta)
        Iq_cnv_ref2 = limit_value*sin(theta)
    else
        Id_cnv_ref2 = Id_cnv_ref
        Iq_cnv_ref2 = Iq_cnv_ref
    end
    return Id_cnv_ref2, Iq_cnv_ref2
end

