function limit_output_current(
    limiter::Nothing,
    Id_cnv_ref::ACCEPTED_REAL_TYPES,
    Iq_cnv_ref::ACCEPTED_REAL_TYPES,
)
    return Id_cnv_ref, Iq_cnv_ref
end

function limit_output_current(
    limiter::PSY.InstantaneousOutputCurrentLimiter,
    Id_cnv_ref::ACCEPTED_REAL_TYPES,
    Iq_cnv_ref::ACCEPTED_REAL_TYPES,
)
    d_lim = PSY.get_Id_max(limiter)
    q_lim = PSY.get_Iq_max(limiter)
    Id_cnv_ref2 = clamp(Id_cnv_ref, -d_lim, d_lim)
    Iq_cnv_ref2 = clamp(Iq_cnv_ref, -q_lim, q_lim)
    return Id_cnv_ref2, Iq_cnv_ref2
end

function limit_output_current(
    limiter::PSY.MagnitudeOutputCurrentLimiter,
    Id_cnv_ref::ACCEPTED_REAL_TYPES,
    Iq_cnv_ref::ACCEPTED_REAL_TYPES,
)
    limit_value = PSY.get_I_max(limiter)
    theta = atan(Iq_cnv_ref, Id_cnv_ref)
    if (Id_cnv_ref^2 + Iq_cnv_ref^2)^(1 / 2) > limit_value
        Id_cnv_ref2 = limit_value * cos(theta)
        Iq_cnv_ref2 = limit_value * sin(theta)
    else
        Id_cnv_ref2 = Id_cnv_ref
        Iq_cnv_ref2 = Iq_cnv_ref
    end
    return Id_cnv_ref2, Iq_cnv_ref2
end

function limit_output_current(
    limiter::PSY.SaturationOutputCurrentLimiter,
    Id_cnv_ref::ACCEPTED_REAL_TYPES,
    Iq_cnv_ref::ACCEPTED_REAL_TYPES,
)
    limit_value = PSY.get_I_max(limiter)
    gain = PSY.get_kw(limiter)
    theta = atan(Iq_cnv_ref, Id_cnv_ref)
    if (Id_cnv_ref^2 + Iq_cnv_ref^2)^(1 / 2) > limit_value
        Id_cnv_ref2 = limit_value * cos(theta)
        Iq_cnv_ref2 = limit_value * sin(theta)
    else
        Id_cnv_ref2 = Id_cnv_ref
        Iq_cnv_ref2 = Iq_cnv_ref
    end
    Del_Vv_d = gain * (Id_cnv_ref - Id_cnv_ref2)
    Del_Vv_q = gain * (Iq_cnv_ref - Iq_cnv_ref2)
    return Id_cnv_ref2, Iq_cnv_ref2, Del_Vv_d, Del_Vv_q
end

function limit_output_current(
    limiter::PSY.HybridOutputCurrentLimiter,
    Id_cnv_ref::ACCEPTED_REAL_TYPES,
    Iq_cnv_ref::ACCEPTED_REAL_TYPES,
    ::ACCEPTED_REAL_TYPES,
)
    limit_value = PSY.get_I_max(limiter)
    real_imped = PSY.get_rv(limiter)
    imag_imped = PSY.get_lv(limiter)
    theta = atan(Iq_cnv_ref, Id_cnv_ref)
    use_clamp = get(PSY.get_ext(limiter), "nondifferentiable", true)
    if use_clamp
        if (Id_cnv_ref^2 + Iq_cnv_ref^2)^(1 / 2) > limit_value
            Id_cnv_ref2 = limit_value * cos(theta)
            Iq_cnv_ref2 = limit_value * sin(theta)
        else
            Id_cnv_ref2 = Id_cnv_ref
            Iq_cnv_ref2 = Iq_cnv_ref
        end
        Del_Vv_d =
            imag_imped * (Iq_cnv_ref2 - Iq_cnv_ref) +
            real_imped * (Id_cnv_ref - Id_cnv_ref2)
        Del_Vv_q =
            imag_imped * (Id_cnv_ref - Id_cnv_ref2) +
            real_imped * (Iq_cnv_ref - Iq_cnv_ref2)
        return Id_cnv_ref2, Iq_cnv_ref2, Del_Vv_d, Del_Vv_q
    else
        ε = 0.01
        Idq = sqrt(Id_cnv_ref^2 + Iq_cnv_ref^2)
        ρ = -ε * log(exp(-1 / ε) + exp(-limit_value / (ε * Idq)))
        Id_cnv_ref2 = ρ * Id_cnv_ref
        Iq_cnv_ref2 = ρ * Iq_cnv_ref
        Del_Vv_d =
            imag_imped * (Iq_cnv_ref2 - Iq_cnv_ref) +
            real_imped * (Id_cnv_ref - Id_cnv_ref2)
        Del_Vv_q =
            imag_imped * (Id_cnv_ref - Id_cnv_ref2) +
            real_imped * (Iq_cnv_ref - Iq_cnv_ref2)
        return Id_cnv_ref2, Iq_cnv_ref2, Del_Vv_d, Del_Vv_q
    end
end
