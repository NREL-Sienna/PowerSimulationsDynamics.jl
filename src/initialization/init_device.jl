function initialize_device(device::PSY.DynamicGenerator)
    #Obtain States
    device_states = zeros(PSY.get_n_states(device))

    #Initialize Machine and Shaft: δ and ω
    initialize_mach_shaft!(device_states, device)
    #Initialize extra Shaft states
    initialize_shaft!(device_states, device)
    #Initialize AVR
    initialize_avr!(device_states, device)
    #Initialize TG
    initialize_tg!(device_states, device)
    #Initialize PSS
    initialize_pss!(device_states, device)

    return device_states
end

function initialize_device(device::PSY.DynamicInverter)
    #Obtain States
    device_states = zeros(PSY.get_n_states(device))

    #Initialize Machine and Shaft: V and I
    initialize_filter!(device_states, device)
    #Initialize freq estimator
    initialize_frequency_estimator!(device_states, device)
    #Initialize OuterLoop
    initialize_outer!(device_states, device)
    #Initialize DCside
    initialize_DCside!(device_states, device)
    #Initialize InnerLoop
    initialize_inner!(device_states, device)
    return device_states
end
