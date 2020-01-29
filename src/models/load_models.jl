function mdl_Zload!(voltage_r,
                    voltage_i,
                    current_r,
                    current_i,
                    device::PSY.PowerLoad,
                    sys::PSY.System)

    #Get parameters
    P = PSY.get_activepower(device)
    Q = PSY.get_reactivepower(device)

    #Compute impedance and load
    #Z_load =  1.0/(P-Q*1im) #in pu
    #I = -(voltage_r[1] + voltage_i[1]*1im)/Z_load #in pu flowing out
    #current_r[1] += real(I)
    #current_i[1] += imag(I)

    #Compute current given constant power
    #I = (P-Q*1im)/(voltage_r[1] - voltage_i[1]*1im)
    #current_r[1] += real(I)
    #current_i[1] += imag(I)

    #Update current
    #This model creates an equivalent RL/RC circuit assuming nominal voltage (1 pu)
    current_r[1] += -(voltage_r[1]*P + voltage_i[1]*Q) #in system pu flowing out
    current_i[1] += -(voltage_i[1]*P - voltage_r[1]*Q) #in system pu flowing out

    return
end
