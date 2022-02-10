"""
Function to obtain the output current time series of a Dynamic Inverter model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific current.

"""
function compute_output_current(
    res::SimulationResults,
    dynamic_device::G,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    dt::Union{Nothing, Float64}
) where {G <: PSY.DynamicInverter}

    #Obtain Data
    sys = get_system(res)

    #Get machine
    filt = PSY.get_filter(dynamic_device)
    converter = PSY.get_converter(dynamic_device)
    Sbase = PSY.get_base_power(sys)
    basepower = PSY.get_base_power(dynamic_device)
    base_power_ratio = basepower / Sbase
    return _output_current(
        filt,
        converter,
        PSY.get_name(dynamic_device),
        V_R,
        V_I,
        base_power_ratio,
        res,
        dynamic_device,
        dt
    )
end

"""
Function to obtain the output current time series of a LCL Filter model out of the DAE Solution. It is dispatched via the Filter type.

"""
function _output_current(
    ::PSY.LCLFilter,
    ::C,
    name::String,
    ::Vector{Float64},
    ::Vector{Float64},
    base_power_ratio::Float64,
    res::SimulationResults,
    dynamic_device::G,
    dt::Union{Nothing, Float64}
) where {C <: PSY.Converter, G <: PSY.DynamicInverter}
    ts, ir_filter = post_proc_state_series(res, (name, :ir_filter), dt)
    ts, ii_filter = post_proc_state_series(res, (name, :ii_filter), dt)

    return ts, base_power_ratio * ir_filter, base_power_ratio * ii_filter
end

"""
Function to obtain the output current time series of a REGCA converter model out of the DAE Solution. It is dispatched via the Converter type.

"""
function _output_current(
    filt::PSY.RLFilter,
    converter::PSY.RenewableEnergyConverterTypeA,
    name::String,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    base_power_ratio::Float64,
    res::SimulationResults,
    ::G,
    dt::Union{Nothing, Float64}
) where {G <: PSY.DynamicInverter}

    #Get Converter parameters
    Lv_pnt0, Lv_pnt1 = PSY.get_Lv_pnts(converter)
    K_hv = PSY.get_K_hv(converter)
    Vo_lim = PSY.get_Vo_lim(converter)
    Brkpt = PSY.get_Brkpt(converter)
    Zerox = PSY.get_Zerox(converter)
    Lvpl1 = PSY.get_Lvpl1(converter)
    Lvpl_sw = PSY.get_Lvpl_sw(converter)
    R_source = PSY.get_R_source(converter)
    X_source = PSY.get_X_source(converter)
    #Get Filter parameters
    rf = PSY.get_rf(filt)
    lf = PSY.get_lf(filt)

    #Compute additional terms
    V_t = sqrt.(V_R .^ 2 + V_I .^ 2)
    G_lv = get_LV_current_gain.(V_t, Lv_pnt0, Lv_pnt1)
    Iq_extra = max.(K_hv * (V_t .- Vo_lim), 0.0)

    #Compute current
    ts, Ip = post_proc_state_series(res, (name, :Ip), dt)
    _, Iq = post_proc_state_series(res, (name, :Iq), dt)
    _, Vmeas = post_proc_state_series(res, (name, :Vmeas), dt)
    Ip_sat = Ip
    if Lvpl_sw == 1
        LVPL = get_LVPL_gain.(Vmeas, Zerox, Brkpt, Lvpl1)
        Ip_sat = Ip <= LVPL ? Ip : LVPL
    end
    Id_cnv = G_lv .* Ip_sat
    Iq_cnv = -Iq - Iq_extra
    θ = atan.(V_I ./ V_R)
    Id_cnv = G_lv .* Ip_sat
    Iq_cnv = -Iq - Iq_extra
    Ir_cnv = Id_cnv .* cos.(θ) - Iq_cnv .* sin.(θ)
    Ii_cnv = Id_cnv .* sin.(θ) + Iq_cnv .* cos.(θ)

    #Compute filter current
    Zf = rf + lf * 1im
    Z_source = R_source + X_source * 1im

    function V_cnv_calc(Ir_cnv, Ii_cnv, Vr_inv, Vi_inv)
        I_cnv = Ir_cnv + 1im * Ii_cnv
        V_inv = Vr_inv + 1im * Vi_inv
        if lf != 0.0 || rf != 0.0
            V_cnv = (I_cnv + V_inv / Zf) / (1.0 / Z_source + 1.0 / Zf)
        else
            V_cnv = V_inv
        end
        return V_cnv
    end

    V_cnv = V_cnv_calc(Ir_cnv, Ii_cnv, V_R, V_I)
    Vr_cnv = real(V_cnv)
    Vi_cnv = imag(V_cnv)
    if lf != 0.0 || rf != 0.0
        Zmag = rf^2 + lf^2
        Ir_filt = (1.0 / Zmag) * ((Vr_cnv - V_R) * rf + (Vi_cnv - V_I) * lf)
        Ii_filt = (1.0 / Zmag) * ((Vi_cnv - V_I) * rf - (Vr_cnv - V_I) * lf)
    else
        Ir_filt = Ir_cnv
        Ii_filt = Ii_cnv
    end

    return ts, base_power_ratio * Ir_filt, base_power_ratio * Ii_filt
end
