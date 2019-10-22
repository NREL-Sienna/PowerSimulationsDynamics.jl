abstract type DCSource <: InverterComponent end

@def dcsource_ports begin
    state_input = Vector{Symbol}()
    inner_input = Vector{Int64}()
end


"""
Parameters of a Fixed DC Source that returns a fixed DC voltage
# Conmutable structor
```julia
FixedDCSource(voltage)
```

# Arguments
* `voltage`::Float64 : Fixed DC voltage to the converter

"""
mutable struct FixedDCSource <: DCSource
    voltage::Float64
    n_states::Int64
    states::Vector{Symbol}
    ports::Ports

        function FixedDCSource(voltage::Float64)

            @dcsource_ports

            new(voltage,
                0,
                Vector{Symbol}(),
                Ports(state_input, inner_input))
        end
end
