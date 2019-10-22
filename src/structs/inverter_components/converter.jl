abstract type Converter <: InverterComponent end

#TODO: converter_ports for non-fixed DC converter will need:
## VDC as state of DC source
## IDC as algebraic inner_var output from cnv real output power


@def converter_ports begin
    state_input = Vector{Symbol}()
    inner_input = [md_var, mq_var, Vdc_var, Vdcnv_var,Vqcnv_var]
end

"""
Parameters of an average converter model
# Conmutable structor
```julia
AvgCnvFixedDC(v_rated, s_rated)
```

# Arguments
*  `v_rated`::Float64 : rated voltage
*  `s_rated`::Float64 : rated VA
"""
mutable struct AvgCnvFixedDC <: Converter
    v_rated::Float64
    s_rated::Float64
    n_states::Int64
    states::Vector{Symbol}
    ports::Ports

        function AvgCnvFixedDC(v_rated::Float64,
                                s_rated::Float64)

            @converter_ports

            new(v_rated,
                s_rated,
                0,
                Vector{Symbol}(),
                Ports(state_input, inner_input))
        end
end
