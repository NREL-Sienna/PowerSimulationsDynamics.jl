function initialize_device!(
    initial_guess,
    ix_range,
    device::PSY.DynamicGenerator
    )
    #Obtain States
    device_states = @view initial_guess[ix_range]
    
    initialize_machine!(device_states, device)
    initialize_shaft!(device_states, device)
    initialize_avr!(device_states, device)
    initialize_tg!(device_states, device)
    initialize_pss!(device_states, device)

end