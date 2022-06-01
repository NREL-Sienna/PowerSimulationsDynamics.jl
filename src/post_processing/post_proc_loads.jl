"""
Function to obtain the output current time series of a 5th Order Induction Machine model out of the DAE Solution. 
It receives the simulation inputs, the dynamic device and bus voltage. 

"""
function compute_output_current(
    res::SimulationResults,
    dynamic_device::PSY.SingleCageInductionMachine,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    dt::Union{Nothing, Float64},
)
    #Obtain Data
    sys = get_system(res)
    name = PSY.get_name(dynamic_device)
    Sbase = PSY.get_base_power(sys)
    basepower = PSY.get_base_power(dynamic_device)
    base_power_ratio = basepower / Sbase

    #Get parameters
    X_ls = PSY.get_X_ls(dynamic_device)
    X_lr = PSY.get_X_lr(dynamic_device)
    X_ad = PSY.get_X_ad(dynamic_device)
    X_aq = PSY.get_X_aq(dynamic_device)
    B_sh = PSY.get_B_shunt(dynamic_device)

    #voltages in Stator
    v_qs = V_I
    v_ds = V_R

    # Read states
    ts, ψ_qs = post_proc_state_series(res, (name, :ψ_qs), dt)
    _, ψ_ds = post_proc_state_series(res, (name, :ψ_ds), dt)
    _, ψ_qr = post_proc_state_series(res, (name, :ψ_qr), dt)
    _, ψ_dr = post_proc_state_series(res, (name, :ψ_dr), dt)

    #Additional Fluxes 
    ψ_mq = X_aq * (ψ_qs / X_ls + ψ_qr / X_lr) # (4.14-15) in Krause
    ψ_md = X_ad * (ψ_ds / X_ls + ψ_dr / X_lr) # (4.14-16) in Krause

    # Stator motor currents in QD 
    i_qs = 1 / X_ls * (ψ_qs - ψ_mq) # (4.14-1) in Krause
    i_ds = 1 / X_ls * (ψ_ds - ψ_md) # (4.14-2) in Krause

    I_R = -base_power_ratio * (i_ds - v_qs * B_sh)
    I_I = -base_power_ratio * (i_qs + v_ds * B_sh)

    return ts, I_R, I_I
end

"""
Function to obtain the output current time series of a 3th Order Induction Machine model out of the DAE Solution. 
It receives the simulation inputs, the dynamic device and bus voltage. 

"""
function compute_output_current(
    res::SimulationResults,
    dynamic_device::PSY.SimplifiedSingleCageInductionMachine,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    dt::Union{Nothing, Float64},
)
    #Obtain Data
    sys = get_system(res)
    name = PSY.get_name(dynamic_device)
    Sbase = PSY.get_base_power(sys)
    basepower = PSY.get_base_power(dynamic_device)
    base_power_ratio = basepower / Sbase

    # Get Parameters
    R_s = PSY.get_R_s(dynamic_device)
    X_m = PSY.get_X_m(dynamic_device)
    B_sh = PSY.get_B_shunt(dynamic_device)
    X_rr = PSY.get_X_rr(dynamic_device)
    X_p = PSY.get_X_p(dynamic_device)

    # Read states
    ts, ψ_qr = post_proc_state_series(res, (name, :ψ_qr), dt)
    _, ψ_dr = post_proc_state_series(res, (name, :ψ_dr), dt)

    # voltages in QD 
    v_qs = V_I
    v_ds = V_R

    # TODO: Read proper sys_ω
    sys_ω = 1.0

    # Stator and rotor currents in QD 
    i_qs =
        1 / (R_s^2 + (sys_ω * X_p)^2) * (
            (R_s * v_qs - sys_ω * X_p * v_ds) -
            (R_s * sys_ω * X_m / X_rr * ψ_dr + sys_ω * X_p * sys_ω * X_m / X_rr * ψ_qr)
        )
    i_ds =
        1 / (R_s^2 + (sys_ω * X_p)^2) * (
            (R_s * v_ds + sys_ω * X_p * v_qs) -
            (-R_s * sys_ω * X_m / X_rr * ψ_qr + sys_ω * X_p * sys_ω * X_m / X_rr * ψ_dr)
        )

    I_R = -base_power_ratio * (i_ds - v_qs * B_sh)
    I_I = -base_power_ratio * (i_qs + v_ds * B_sh)

    return ts, I_R, I_I
end
