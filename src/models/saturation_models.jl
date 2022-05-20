"""
    Saturation function for quadratic saturation models for machines
        Se(x) = B * (x - A)^2 / x
"""
function saturation_function(
    machine::Union{PSY.RoundRotorQuadratic, PSY.SalientPoleQuadratic},
    x::ACCEPTED_REAL_TYPES,
)
    Sat_A, Sat_B = PSY.get_saturation_coeffs(machine)
    return Sat_B * (x - Sat_A)^2 / x
end

"""
    Saturation function for exponential saturation models for machines
        Se(x) = B * x^A
"""
function saturation_function(
    machine::Union{PSY.RoundRotorExponential, PSY.SalientPoleExponential},
    x::ACCEPTED_REAL_TYPES,
)
    Sat_A, Sat_B = PSY.get_saturation_coeffs(machine)
    return Sat_B * x^Sat_A
end

function rectifier_function(I::Float64)
    if I <= 0.0
        return 1.0
    elseif I <= 0.433
        return 1.0 - 0.577 * I
    elseif I < 0.75
        return sqrt(0.75 - I^2)
    elseif I <= 1.0
        return 1.732 * (1.0 - I)
    else
        return 0.0
    end
end

function saturation_function(avr::Union{PSY.ESAC1A, PSY.EXAC1}, x::ACCEPTED_REAL_TYPES)
    Sat_A, Sat_B = PSY.get_saturation_coeffs(avr)
    return Sat_B * (x - Sat_A)^2 / x
end

function rectifier_function(I::T) where {T <: ACCEPTED_REAL_TYPES}
    if I <= 0.0
        return one(T)
    elseif I <= 0.433
        return 1.0 - 0.577 * I
    elseif I < 0.75
        return sqrt(0.75 - I^2)
    elseif I <= 1.0
        return 1.732 * (1.0 - I)
    else
        return zero(T)
    end
end

function output_pss_limiter(
    V_ss::X,
    V_ct::X,
    V_cl::Float64,
    V_cu::Float64,
) where {X <: ACCEPTED_REAL_TYPES}
    # Bypass limiter block if one limiter parameter is set to zero.
    if V_cl == 0.0 || V_cu == 0.0
        return V_ss
    end
    if V_cl <= V_ct <= V_cu
        return V_ss
    elseif V_ct < V_cl
        return zero(X)
    elseif V_ct > V_cu
        return zero(X)
    else
        error("Logic error in PSS")
    end
end

function deadband_function(
    x::T,
    db_low::Float64,
    db_high::Float64,
) where {T <: ACCEPTED_REAL_TYPES}
    if x > db_high
        return x - db_high
    elseif x < db_low
        return x - db_low
    else
        return zero(T)
    end
end

function current_limit_logic(
    inner_control::PSY.RECurrentControlB,
    ::Val{0}, #PQ_Flag = 0: Q Priority
    Vt_filt::X,
    Ip_cmd::X,
    Iq_cmd::X,
) where {X <: ACCEPTED_REAL_TYPES}
    I_max = PSY.get_I_max(inner_control)
    Iq_max = I_max
    Iq_min = -Iq_max
    Ip_min = 0.0
    local_I = I_max^2 - Iq_cmd^2
    if local_I < 0
        local_I = 0
    else
        local_I = sqrt(local_I)
    end
    if local_I < I_max
        Ip_max = local_I
    else
        Ip_max = I_max
    end
    return Ip_min, Ip_max, Iq_min, Iq_max
end

function current_limit_logic(
    inner_control::PSY.RECurrentControlB,
    ::Val{1}, #PQ_Flag = 1: P Priority
    Vt_filt::X,
    Ip_cmd::X,
    Iq_cmd::X,
) where {X <: ACCEPTED_REAL_TYPES}
    I_max = PSY.get_I_max(inner_control)
    Ip_max = I_max
    Ip_min = 0.0
    local_I = I_max^2 - Ip_cmd^2
    if local_I < 0
        local_I = 0
    else
        local_I = sqrt(local_I)
    end
    if local_I < I_max #fixed typo 
        Iq_max = local_I
    else
        Iq_max = I_max
    end
    Iq_min = -Iq_max
    return Ip_min, Ip_max, Iq_min, Iq_max
end

function current_limit_logic(
    device::PSY.AggregateDistributedGenerationA,
    Ip_cmd::X,
    Iq_cmd::X,
) where {X <: ACCEPTED_REAL_TYPES}
    PQ_Flag = PSY.get_PQ_Flag(device)
    Gen_Flag = PSY.get_Gen_Flag(device)
    I_max = PSY.get_I_max(device)

    if PQ_Flag == 1  # P Priority 
        Ip_max = I_max
        local_I = I_max^2 - Ip_cmd^2
        if local_I < 0
            local_I = 0
        else
            local_I = sqrt(local_I)
        end
        if local_I < I_max
            Iq_max = local_I
        else
            Iq_max = I_max
        end
    elseif PQ_Flag == 0     #Q Priority  
        Iq_max = I_max
        local_I = I_max^2 - Iq_cmd^2
        if local_I < 0
            local_I = 0
        else
            local_I = sqrt(local_I)
        end
        if local_I < I_max
            Ip_max = local_I
        else
            Ip_max = I_max
        end
    else
        @error "Unsupported value of PQ_Flag"
    end
    if Gen_Flag == 1
        Ip_min = 0
    elseif Gen_Flag == 0
        Ip_min = -Ip_max
    else
        @error "Unsupported value of Gen_Flag"
    end
    Iq_min = -Iq_max
    return Ip_min, Ip_max, Iq_min, Iq_max
end

function voltage_trip_logic!(
    inner_vars::AbstractArray{X},
    device::PSY.AggregateDistributedGenerationA,
    Vmeas::X,
    t,
) where {X <: ACCEPTED_REAL_TYPES}
    #Read parameters
    Vtrip_Flag = PSY.get_Vtrip_Flag(device)
    (tvl0, vl0), (tvl1, vl1) = PSY.get_vl_pnts(device)
    (tvh0, vh0), (tvh1, vh1) = PSY.get_vh_pnts(device)
    Vrfrac = PSY.get_Vrfrac(device)

    if Vtrip_Flag == 0
        return 1.0
    end

    #Read inner vars 
    Vmin = inner_vars[Vmin_var]
    Vmax = inner_vars[Vmax_var]
    TimeBelowVl1 = inner_vars[TimeBelowVl1_var]
    TimeAboveVh1 = inner_vars[TimeAboveVh1_var]
    TimeBelowVl0 = inner_vars[TimeBelowVl0_var]
    TimeAboveVh0 = inner_vars[TimeAboveVh0_var]
    ActiveTimerVl1 = inner_vars[ActiveTimerVl1_var]
    ActiveTimerVh1 = inner_vars[ActiveTimerVh1_var]
    ActiveFracLow = inner_vars[ActiveFracLow_var]
    ActiveFracHigh = inner_vars[ActiveFracHigh_var]
    ActiveTimerVl0 = inner_vars[ActiveTimerVl0_var]
    ActiveTimerVh0 = inner_vars[ActiveTimerVh0_var]
    ActiveTripLow = inner_vars[ActiveTripLow_var]
    ActiveTripHigh = inner_vars[ActiveTripHigh_var]

    #TODO - should this happen in init_device and be stored in ext, or OK to implement here?
    if t == 0
        Vmin = vl1
        Vmax = vh1
    end

    #Starting Timers for Low Voltage 
    if Vmeas >= vl1
        ActiveTimerVl1 = 0
    elseif ActiveTimerVl1 == 0
        ActiveTimerVl1 = 1
        TimeBelowVl1 = t
    end
    if Vmeas >= vl0
        ActiveTimerVl0 = 0
    elseif ActiveTimerVl0 == 0
        ActiveTimerVl0 = 1
        TimeBelowVl0 = t
    end

    #Starting Timers for High Voltage
    if Vmeas <= vh1
        ActiveTimerVh1 = 0
    elseif ActiveTimerVh1 == 0
        ActiveTimerVh1 = 1
        TimeAboveVh1 = t
    end
    if Vmeas <= vh0
        ActiveTimerVh0 = 0
    elseif ActiveTimerVh0 == 0
        ActiveTimerVh0 = 1
        TimeAboveVh0 = t
    end

    #Use timers to see if on fractional restart curve 
    if (ActiveFracLow == 0) && (ActiveTimerVl1 == 1) && ((t - TimeBelowVl1) >= tvl1)
        ActiveFracLow = 1
    end
    if (ActiveFracHigh == 0) && (ActiveTimerVh1 == 1) && ((t - TimeAboveVh1) >= tvh1)
        ActiveFracHigh = 1
    end

    #Use timers to see if completely tripped
    if (ActiveTripLow == 0) && (ActiveTimerVl0 == 1) && ((t - TimeBelowVl0) >= tvl0)
        ActiveTripLow = 1
    end
    if (ActiveTripHigh == 0) && (ActiveTimerVh0 == 1) && ((t - TimeAboveVh0) >= tvh0)
        ActiveTripHigh = 1
    end

    #Keep track of the minimum and maximum voltages
    if (Vmin > Vmeas) && (ActiveFracLow == 1)
        Vmin = Vmeas
        if Vmin < vl0
            Vmin = vl0
        end
    end
    if (Vmax < Vmeas) && (ActiveFracHigh == 1)
        Vmax = Vmeas
        if Vmax > vh0
            Vmax = vh0
        end
    end

    #Calculate the low voltage multiplier
    if (Vmeas <= vl0) || (ActiveTripLow == 1)
        Vlmult = 0.0
    elseif (Vmeas <= vl1) && (Vmeas > Vmin) && (ActiveFracLow == 1)
        Vlmult = ((Vmin - vl0) + Vrfrac * (Vmeas - Vmin)) / (vl1 - vl0)
    elseif (Vmeas <= vl1)
        Vlmult = (Vmeas - vl0) / (vl1 - vl0)
    elseif (ActiveFracLow == 0)
        Vlmult = 1.0
    else
        Vlmult = ((Vmin - vl0) + Vrfrac * (vl1 - Vmin)) / (vl1 - vl0)
    end
    #Calculate the high voltage multiplier 
    if (Vmeas >= vh0) || (ActiveTripHigh == 1)
        Vhmult = 0.0
    elseif (Vmeas >= vh1) && (Vmeas < Vmax) && (ActiveFracHigh == 1)
        Vhmult = ((Vmax - vh0) + Vrfrac * (Vmeas - Vmax)) / (vh1 - vh0)
    elseif (Vmeas >= vh1)
        Vhmult = (Vmeas - vh0) / (vh1 - vh0)
    elseif (ActiveFracHigh == 0)
        Vhmult = 1.0
    else
        Vhmult = ((Vmax - vh0) + Vrfrac * (vh1 - Vmax)) / (vh1 - vh0)
    end

    Mult = Vlmult * Vhmult
    #Update all the inner vars 
    inner_vars[Vmin_var] = Vmin
    inner_vars[Vmax_var] = Vmax
    inner_vars[TimeBelowVl1_var] = TimeBelowVl1
    inner_vars[TimeAboveVh1_var] = TimeAboveVh1
    inner_vars[TimeBelowVl0_var] = TimeBelowVl0
    inner_vars[TimeAboveVh0_var] = TimeAboveVh0
    inner_vars[ActiveTimerVl1_var] = ActiveTimerVl1
    inner_vars[ActiveTimerVh1_var] = ActiveTimerVh1
    inner_vars[ActiveFracLow_var] = ActiveFracLow
    inner_vars[ActiveFracHigh_var] = ActiveFracHigh
    inner_vars[ActiveTimerVl0_var] = ActiveTimerVl0
    inner_vars[ActiveTimerVh0_var] = ActiveTimerVh0
    inner_vars[ActiveTripLow_var] = ActiveTripLow
    inner_vars[ActiveTripHigh_var] = ActiveTripHigh

    return Mult
end

function frequency_trip_logic!(
    inner_vars::AbstractArray{X},
    device::PSY.AggregateDistributedGenerationA,
    Fmeas::X,
    Vt::X,
    t,
) where {X <: ACCEPTED_REAL_TYPES}

    #Read parameters
    fl = PSY.get_fl(device)
    fh = PSY.get_fh(device)
    tfl = PSY.get_tfl(device)
    tfh = PSY.get_tfh(device)
    Vpr = PSY.get_Vpr(device)
    Ftrip_Flag = PSY.get_Ftrip_Flag(device)

    if Ftrip_Flag == 0
        return 1.0
    end

    #Read inner vars 
    TimeBelowFl = inner_vars[TimeBelowFl_var]
    TimeAboveFh = inner_vars[TimeAboveFh_var]
    ActiveTimerFl = inner_vars[ActiveTimerFl_var]
    ActiveTimerFh = inner_vars[ActiveTimerFh_var]
    ActiveTripFLow = inner_vars[ActiveTripFLow_var]
    ActiveTripFHigh = inner_vars[ActiveTripFHigh_var]

    #Calculate Frelay 
    if Vt > Vpr
        Frelay = Fmeas
    else
         Frelay = 1.0 
    end 

    #Start timers 
    if Frelay >= fl
        ActiveTimerFl = 0
    elseif ActiveTimerFl == 0
        ActiveTimerFl = 1
        TimeBelowFl = t
    end
    if Frelay <= fh
        ActiveTimerFh = 0
    elseif ActiveTimerFh == 0
        ActiveTimerFh = 1
        TimeAboveFh = t
    end

    #Use timers to see if tripped 
    if (ActiveTripFLow == 0) && (ActiveTimerFl == 1) && ((t - TimeBelowFl) >= tfl)
        ActiveTripFLow = 1
    end
    if (ActiveTripFHigh == 0) && (ActiveTimerFh == 1) && ((t - TimeAboveFh) >= tfh)
        ActiveTripFHigh = 1
    end

    #Calculate Fmult 
    if (ActiveTripFLow == 1) || (ActiveTripFHigh == 1)
        Fmult = 0.0 
    else 
        Fmult = 1.0 
    end  

    #Update inner vars 
    inner_vars[TimeBelowFl_var] = TimeBelowFl
    inner_vars[TimeAboveFh_var] = TimeAboveFh
    inner_vars[ActiveTimerFl_var] = ActiveTimerFl
    inner_vars[ActiveTimerFh_var] = ActiveTimerFh
    inner_vars[ActiveTripFLow_var] = ActiveTripFLow
    inner_vars[ActiveTripFHigh_var] = ActiveTripFHigh

    return Fmult
end

function get_LVPL_gain(
    Vmeas::T,
    Zerox::Float64,
    Brkpt::Float64,
    Lvpl1::Float64,
) where {T <: ACCEPTED_REAL_TYPES}
    if Vmeas < Zerox
        return zero(T)
    elseif Vmeas > Brkpt
        return one(T) * Lvpl1
    else
        return Lvpl1 / (Brkpt - Zerox) * (Vmeas - Zerox)
    end
end

function get_LV_current_gain(
    V_t::T,
    Lv_pnt0::Float64,
    Lv_pnt1::Float64,
) where {T <: ACCEPTED_REAL_TYPES}
    if V_t < Lv_pnt0
        return zero(T)
    elseif V_t > Lv_pnt1
        return one(T)
    else
        return 1.0 / (Lv_pnt0 - Lv_pnt1) * (V_t - Lv_pnt0)
    end
end
