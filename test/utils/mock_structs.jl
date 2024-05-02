struct MockIntegrator
    du::Vector
    u::Vector
    p::Any
end

function MockIntegrator(inputs)
    return MockIntegrator(
        zeros(PSID.get_variable_count(inputs)),
        zeros(PSID.get_variable_count(inputs)),
        PSID.get_parameters(inputs),
    )
end
