function mdl_source!(voltage_r,
                     voltage_i,
                        current_r,
                        current_i,
                     device::StaticSource,
                     sys::DynamicSystem)


        #Load device parameters
        V_R = device.V_R
        V_I = device.V_I
        X_th = device.X_th

        #I = ( (V_R + V_I*1im) - (V_tR + V_tI*1im) )/(X_th*1im)
        #I_r = (V_I - V_tI)/X_th #in system pu flowing out
        #I_i = -(V_R - V_tR)/X_th #in system pu flowing out

        #update current
        current_r[1] += (V_I - voltage_i[1])/X_th #in system pu flowing out
        current_i[1] += -(V_R - voltage_r[1])/X_th #in system pu flowing out

        return
end
