abstract type Filter <: InverterComponent end

@def filter_ports begin
    #TODO: If converter has dynamics, need to connect state_input
    state_input = [:δθ_vsm] #[:Vd_c, :Vq_c] #, :Id_c, :Iq_c]
    inner_input = [VR_inv_var,VI_inv_var,Vdcnv_var,Vqcnv_var,δdqRI_var,Vdo_var,Vqo_var]
end

"""
*  `lf`::Float64 : filter inductance
*  `rf`::Float64 : filter resistance
*  `cf`::Float64 : filter capacitance
"""
struct LCFilter <: Filter
    lf::Float64
    rf::Float64
    cf::Float64
    n_states::Int64
    states::Vector{Symbol}
    ports::Ports

        function LCFilter(lf::Float64,
                          rf::Float64,
                          cf::Float64)

            @filter_ports

            new(lf,
                rf,
                cf,
                4,
                [:id_o, :iq_o],
                Ports(state_input, inner_input))

        end

    end


"""
Parameters of a LCL filter outside the converter
# Conmutable structor
```julia
LCLFilter(lf, rf, cf, lg, rg)
```

# Arguments
* `lf`::Float64 : Series inductance in p.u. of converter filter
* `rf`::Float64 : Series resistance in p.u. of converter filter
* `cf`::Float64 : Shunt capacitance in p.u. of converter filter
* `lg`::Float64 : Series inductance in p.u. of converter filter to the grid
* `rg`::Float64 : Series resistance in p.u. of converter filter to the grid

"""
struct LCLFilter <: Filter
        lf::Float64
        rf::Float64
        cf::Float64
        lg::Float64
        rg::Float64
        n_states::Int64
        states::Vector{Symbol}
        ports::Ports

            function LCLFilter(lf::Float64,
                               rf::Float64,
                               cf::Float64,
                               lg::Float64,
                               rg::Float64)

                n_states = 6
                states = [:id_c, :iq_c, :vd_cap, :vq_cap, :id_o, :iq_o]
                @filter_ports

                new(lf,
                    rf,
                    cf,
                    lg,
                    rg,
                    n_states,
                    states,
                    Ports(state_input, inner_input))
            end
end
