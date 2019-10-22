abstract type Machine <: GeneratorComponent end

@def machine_ports begin
    state_input = [:δ, :ω, :Vf]
    inner_input = [VI_gen_var, VR_gen_var]
end


"""
Parameters of 0-states synchronous machine: Classical Model

# Conmutable structor
```julia
BaseMachine(R, Xd_p, eq_p, MVABase)
```

# Arguments
* `R`::Float64 : Resistance after EMF in machine per unit
* `Xd_p`::Float64 : Reactance after EMF in machine per unit
* `eq_p`::Float64 : Fixed EMF behind the impedance
* `MVABase`::Float64 : Nominal Capacity in MVA

"""
mutable struct BaseMachine <: Machine
  R::Float64
  Xd_p::Float64
  eq_p::Float64
  MVABase::Float64
  n_states::Int64
  states::Vector{Symbol}
  ports::Ports

    function BaseMachine(R::Float64,
                             Xd_p::Float64,
                             eq_p::Float64,
                             MVABase::Float64)
            n_states = 0
            states = Vector{Symbol}()
            @machine_ports

        new(R,
            Xd_p,
            eq_p,
            MVABase,
            n_states,
            states,
            Ports(state_input, inner_input))
    end
end


"""
Parameters of 2-states synchronous machine: One d- and One q-Axis Model
# Conmutable structor
```julia
OneDOneQMachine(R, Xd, Xq, Xd_p, Xq_pp, Td0_p, Tq0_pp, MVABase)
```

# Arguments
* `R`::Float64
* `Xd`::Float64
* `Xq`::Float64
* `Xd_p`::Float64 : Reactance after EMF in machine per unit
* `Xq_pp`::Float64
* `Td0_p`::Float64 : Time constant of transient d-axis voltage
* `Tq0_pp`::Float64 : Time constant of subtransient q-axis voltage
* `MVABase`::Float64 : Nominal Capacity in MVA

"""
mutable struct OneDOneQMachine <: Machine
  R::Float64
  Xd::Float64
  Xq::Float64
  Xd_p::Float64
  Xq_p::Float64
  Td0_p::Float64
  Tq0_p::Float64
  MVABase::Float64
  n_states::Int64
  states::Vector{Symbol}
  ports::Ports


      function OneDOneQMachine(R::Float64,
                                Xd::Float64,
                                Xq::Float64,
                                Xd_p::Float64,
                                Xq_p::Float64,
                                Td0_p::Float64,
                                Tq0_p::Float64,
                                MVABase::Float64)
            n_states = 2
            states = [:eq_p, :ed_p]
            @machine_ports

        new(R,
            Xd,
            Xq,
            Xd_p,
            Xq_p,
            Td0_p,
            Tq0_p,
            MVABase,
            n_states,
            states,
            Ports(state_input, inner_input))
    end
end


"""
Parameters of 6-states synchronous machine: Marconato model
# Conmutable structor
```julia
MarconatoMachine(R, Xd, Xq, Xd_p, Xq_p, Xd_pp, Xq_pp, Td0_p, Tq0_p, Td0_pp, Tq0_pp, T_AA, MVABase)
```

# Arguments
* `R`::Float64 : Resistance after EMF in machine per unit
* `Xd`::Float64 : Reactance after EMF in d-axis per unit
* `Xq`::Float64 : Reactance after EMF in q-axis per unit
* `Xd_p`::Float64 : Transient reactance after EMF in d-axis per unit
* `Xq_p`::Float64 : Transient reactance after EMF in q-axis per unit
* `Xd_pp`::Float64 : Subtransient reactance after EMF in d-axis per unit
* `Xq_pp`::Float64 : Subtransient reactance after EMF in q-axis per unit
* `Td0_p`::Float64 : Time constant of transient d-axis voltage
* `Tq0_p`::Float64 : Time constant of transient q-axis voltage
* `Td0_pp`::Float64 : Time constant of subtransient d-axis voltage
* `Tq0_pp`::Float64 : Time constant of subtransient q-axis voltage
* `T_AA`::Float64 : Time constant of d-axis additional leakage
* `MVABase`::Float64 : Nominal Capacity in MVA

"""
mutable struct MarconatoMachine <: Machine
  R::Float64
  Xd::Float64
  Xq::Float64
  Xd_p::Float64
  Xq_p::Float64
  Xd_pp::Float64
  Xq_pp::Float64
  Td0_p::Float64
  Tq0_p::Float64
  Td0_pp::Float64
  Tq0_pp::Float64
  T_AA::Float64
  γd::Float64
  γq::Float64
  MVABase::Float64
  n_states::Int64
  states::Vector{Symbol}
  ports::Ports

      function MarconatoMachine(R::Float64,
                                Xd::Float64,
                                Xq::Float64,
                                Xd_p::Float64,
                                Xq_p::Float64,
                                Xd_pp::Float64,
                                Xq_pp::Float64,
                                Td0_p::Float64,
                                Tq0_p::Float64,
                                Td0_pp::Float64,
                                Tq0_pp::Float64,
                                T_AA::Float64,
                                MVABase::Float64)
            n_states = 6
            states = [:ψq, :ψd, :eq_p, :ed_p, :eq_pp, :ed_pp]
            γd = ( (Td0_pp*Xd_pp)/(Td0_p*Xd_p) )*(Xd-Xd_p)
            γq = ( (Tq0_pp*Xq_pp)/(Tq0_p*Xq_p) )*(Xq-Xq_p)
            @machine_ports

        new(R,
            Xd,
            Xq,
            Xd_p,
            Xq_p,
            Xd_pp,
            Xq_pp,
            Td0_p,
            Tq0_p,
            Td0_pp,
            Tq0_pp,
            T_AA,
            γd,
            γq,
            MVABase,
            n_states,
            states,
            Ports(state_input, inner_input))
    end
end


"""
Parameters of 4-states synchronous machine: Simplified Marconato model
The derivative of stator fluxes (ψd and ψq) is neglected and ωψd = ψd and
ωψq = ψq is assumed (i.e. ω=1.0). This is standard when transmission network
dynamics is neglected.
# Conmutable structor
```julia
SimpleMarconatoMachine(R, Xd, Xq, Xd_p, Xq_p, Xd_pp, Xq_pp, Td0_p, Tq0_p, Td0_pp, Tq0_pp, T_AA, MVABase)
```

# Arguments
* `R`::Float64 : Resistance after EMF in machine per unit
* `Xd`::Float64 : Reactance after EMF in d-axis per unit
* `Xq`::Float64 : Reactance after EMF in q-axis per unit
* `Xd_p`::Float64 : Transient reactance after EMF in d-axis per unit
* `Xq_p`::Float64 : Transient reactance after EMF in q-axis per unit
* `Xd_pp`::Float64 : Subtransient reactance after EMF in d-axis per unit
* `Xq_pp`::Float64 : Subtransient reactance after EMF in q-axis per unit
* `Td0_p`::Float64 : Time constant of transient d-axis voltage
* `Tq0_p`::Float64 : Time constant of transient q-axis voltage
* `Td0_pp`::Float64 : Time constant of subtransient d-axis voltage
* `Tq0_pp`::Float64 : Time constant of subtransient q-axis voltage
* `T_AA`::Float64 : Time constant of d-axis additional leakage
* `MVABase`::Float64 : Nominal Capacity in MVA

"""
mutable struct SimpleMarconatoMachine <: Machine
  R::Float64
  Xd::Float64
  Xq::Float64
  Xd_p::Float64
  Xq_p::Float64
  Xd_pp::Float64
  Xq_pp::Float64
  Td0_p::Float64
  Tq0_p::Float64
  Td0_pp::Float64
  Tq0_pp::Float64
  T_AA::Float64
  γd::Float64
  γq::Float64
  MVABase::Float64
  n_states::Int64
  states::Vector{Symbol}
  ports::Ports

      function SimpleMarconatoMachine(R::Float64,
                                      Xd::Float64,
                                      Xq::Float64,
                                      Xd_p::Float64,
                                      Xq_p::Float64,
                                      Xd_pp::Float64,
                                      Xq_pp::Float64,
                                      Td0_p::Float64,
                                      Tq0_p::Float64,
                                      Td0_pp::Float64,
                                      Tq0_pp::Float64,
                                      T_AA::Float64,
                                      MVABase::Float64)
            n_states = 4
            states = [:eq_p, :ed_p, :eq_pp, :ed_pp]
            γd = ( (Td0_pp*Xd_pp)/(Td0_p*Xd_p) )*(Xd-Xd_p)
            γq = ( (Tq0_pp*Xq_pp)/(Tq0_p*Xq_p) )*(Xq-Xq_p)
            @machine_ports

        new(R,
            Xd,
            Xq,
            Xd_p,
            Xq_p,
            Xd_pp,
            Xq_pp,
            Td0_p,
            Tq0_p,
            Td0_pp,
            Tq0_pp,
            T_AA,
            γd,
            γq,
            MVABase,
            n_states,
            states,
            Ports(state_input, inner_input))
    end
end



"""
Parameters of 6-states synchronous machine: Anderson-Fouad model
# Conmutable structor
```julia
AndersonFouadMachine(R, Xd, Xq, Xd_p, Xq_p, Xd_pp, Xq_pp, Td0_p, Tq0_p, Td0_pp, Tq0_pp, MVABase)
```

# Arguments
* `R`::Float64 : Resistance after EMF in machine per unit
* `Xd`::Float64 : Reactance after EMF in d-axis per unit
* `Xq`::Float64 : Reactance after EMF in q-axis per unit
* `Xd_p`::Float64 : Transient reactance after EMF in d-axis per unit
* `Xq_p`::Float64 : Transient reactance after EMF in q-axis per unit
* `Xd_pp`::Float64 : Subtransient reactance after EMF in d-axis per unit
* `Xq_pp`::Float64 : Subtransient reactance after EMF in q-axis per unit
* `Td0_p`::Float64 : Time constant of transient d-axis voltage
* `Tq0_p`::Float64 : Time constant of transient q-axis voltage
* `Td0_pp`::Float64 : Time constant of subtransient d-axis voltage
* `Tq0_pp`::Float64 : Time constant of subtransient q-axis voltage
* `MVABase`::Float64 : Nominal Capacity in MVA

"""
mutable struct AndersonFouadMachine <: Machine
  R::Float64
  Xd::Float64
  Xq::Float64
  Xd_p::Float64
  Xq_p::Float64
  Xd_pp::Float64
  Xq_pp::Float64
  Td0_p::Float64
  Tq0_p::Float64
  Td0_pp::Float64
  Tq0_pp::Float64
  MVABase::Float64
  n_states::Int64
  states::Vector{Symbol}
  ports::Ports

      function AndersonFouadMachine(R::Float64,
                                Xd::Float64,
                                Xq::Float64,
                                Xd_p::Float64,
                                Xq_p::Float64,
                                Xd_pp::Float64,
                                Xq_pp::Float64,
                                Td0_p::Float64,
                                Tq0_p::Float64,
                                Td0_pp::Float64,
                                Tq0_pp::Float64,
                                MVABase::Float64)
            n_states = 6
            states = [:ψq, :ψd, :eq_p, :ed_p, :eq_pp, :ed_pp]
            @machine_ports

        new(R,
            Xd,
            Xq,
            Xd_p,
            Xq_p,
            Xd_pp,
            Xq_pp,
            Td0_p,
            Tq0_p,
            Td0_pp,
            Tq0_pp,
            MVABase,
            n_states,
            states,
            Ports(state_input, inner_input))
    end
end


"""

Parameters of 4-states simplified Anderson-Fouad (SimpleAFMachine) model.
The derivative of stator fluxes (ψd and ψq) is neglected and ωψd = ψd and
ωψq = ψq is assumed (i.e. ω=1.0). This is standard when transmission network
dynamics is neglected.
If transmission dynamics is considered use the full order Anderson Fouad model.

# Conmutable structor
```julia
SimpleAFMachine(R, Xd, Xq, Xd_p, Xq_p, Xd_pp, Xq_pp, Td0_p, Tq0_p, Td0_pp, Tq0_pp, MVABase)
```

# Arguments
* `R`::Float64 : Resistance after EMF in machine per unit
* `Xd`::Float64 : Reactance after EMF in d-axis per unit
* `Xq`::Float64 : Reactance after EMF in q-axis per unit
* `Xd_p`::Float64 : Transient reactance after EMF in d-axis per unit
* `Xq_p`::Float64 : Transient reactance after EMF in q-axis per unit
* `Xd_pp`::Float64 : Subtransient reactance after EMF in d-axis per unit
* `Xq_pp`::Float64 : Subtransient reactance after EMF in q-axis per unit
* `Td0_p`::Float64 : Time constant of transient d-axis voltage
* `Tq0_p`::Float64 : Time constant of transient q-axis voltage
* `Td0_pp`::Float64 : Time constant of subtransient d-axis voltage
* `Tq0_pp`::Float64 : Time constant of subtransient q-axis voltage
* `MVABase`::Float64 : Nominal Capacity in MVA

"""
mutable struct SimpleAFMachine <: Machine
  R::Float64
  Xd::Float64
  Xq::Float64
  Xd_p::Float64
  Xq_p::Float64
  Xd_pp::Float64
  Xq_pp::Float64
  Td0_p::Float64
  Tq0_p::Float64
  Td0_pp::Float64
  Tq0_pp::Float64
  MVABase::Float64
  n_states::Int64
  states::Vector{Symbol}
  ports::Ports

      function SimpleAFMachine(R::Float64,
                                Xd::Float64,
                                Xq::Float64,
                                Xd_p::Float64,
                                Xq_p::Float64,
                                Xd_pp::Float64,
                                Xq_pp::Float64,
                                Td0_p::Float64,
                                Tq0_p::Float64,
                                Td0_pp::Float64,
                                Tq0_pp::Float64,
                                MVABase::Float64)
            n_states = 4
            states = [:eq_p, :ed_p, :eq_pp, :ed_pp]
            @machine_ports

        new(R,
            Xd,
            Xq,
            Xd_p,
            Xq_p,
            Xd_pp,
            Xq_pp,
            Td0_p,
            Tq0_p,
            Td0_pp,
            Tq0_pp,
            MVABase,
            n_states,
            states,
            Ports(state_input, inner_input))
    end
end




"""

Parameter of a full order flux stator-rotor model without zero sequence flux in the stator.
The derivative of stator fluxes (ψd and ψq) is NOT neglected. Only one q-axis damping circuit
is considered. All parameters are in machine per unit.

Refer to Chapter 3 of Power System Stability and Control by P. Kundur or Chapter 11 of Power
System Dynamics: Stability and Control, by J. Machowski, J. Bialek and J. Bumby, for more details.
Note that the models are somewhat different (but equivalent) due to the different Park
Transformation used in both books.

# Conmutable structor
```julia
FullMachine(R, R_f, R_1d, R_1q, L_d, L_q, L_ad, L_aq, L_f1d, L_f, L_1d, L_1q)
```

# Arguments
* `R`::Float64 : Stator resistance after EMF in per unit
* `R_f`::Float64 : Field rotor winding resistance in per unit
* `R_1d`::Float64 : Damping rotor winding resistance on d-axis in per unit.
                  This value is denoted as RD in Machowski.
* `R_1q`::Float64 : Damping rotor winding resistance on q-axis in per unit.
                  This value is denoted as RQ in Machowski.
* `L_d`::Float64 : Inductance of fictitious damping that represent the effect
                  of the three-phase stator winding in the d-axis of the rotor, in per unit.
                  This value is denoted as L_ad + L_l in Kundur (and Ld in Machowski).
* `L_q`::Float64 : Inductance of fictitious damping that represent the effect
                  of the three-phase stator winding in the q-axis of the rotor, in per unit.
                  This value is denoted as L_aq + L_l in Kundur.
* `L_ad`::Float64 : Mutual inductance between stator winding and rotor field (and damping)
                    winding inductance on d-axis, in per unit
* `L_aq`::Float64 : Mutual inductance between stator winding and rotor damping winding
                    inductance on q-axis, in per unit
* `L_f1d`::Float64 : Mutual inductance between rotor field winding and rotor damping winding
                    inductance on d-axis, in per unit
* `L_ff`::Float64 : Field rotor winding inductance, in per unit
* `L_1d`::Float64 : Inductance of the d-axis rotor damping circuit, in per unit
* `L_1q`::Float64 : Inductance of the q-axis rotor damping circuit, in per unit
* `MVABase`::Float64 : Nominal Capacity in MVA

"""
mutable struct FullMachine <: Machine
  R::Float64
  R_f::Float64
  R_1d::Float64
  R_1q::Float64
  L_d::Float64
  L_q::Float64
  L_ad::Float64
  L_aq::Float64
  L_f1d::Float64
  L_ff::Float64
  L_1d::Float64
  L_1q::Float64
  MVABase::Float64
  inv_d_fluxlink::Array{Float64, 2}
  inv_q_fluxlink::Array{Float64, 2}
  n_states::Int64
  states::Vector{Symbol}
  ports::Ports

      function FullMachine(R::Float64,
                                R_f::Float64,
                                R_1d::Float64,
                                R_1q::Float64,
                                L_d::Float64,
                                L_q::Float64,
                                L_ad::Float64,
                                L_aq::Float64,
                                L_f1d::Float64,
                                L_ff::Float64,
                                L_1d::Float64,
                                L_1q::Float64,
                                MVABase::Float64)

            inv_d_fluxlink = inv( [ [-L_d L_ad L_ad]; #ψd: Eq 3.127 in Kundur
                                    [-L_ad L_ff L_f1d ]; #ψf: Eq 3.130 in Kundur
                                    [-L_ad L_f1d L_1d ] ]) #ψ1d: Eq 3.131 in Kundur
            inv_q_fluxlink = inv( [ [-L_q L_aq]; #ψq: Eq 3.128 in Kundur
                                    [-L_aq L_1q ] ]) #ψ1q: Eq 3.132 in Kundur
            n_states = 5
            states = [:ψd, :ψq, :ψf, :ψ1d, :ψ1q]
            @machine_ports

        new(R,
            R_f,
            R_1d,
            R_1q,
            L_d,
            L_q,
            L_ad,
            L_aq,
            L_f1d,
            L_ff,
            L_1d,
            L_1q,
            MVABase,
            inv_d_fluxlink,
            inv_q_fluxlink,
            n_states,
            states,
            Ports(state_input, inner_input))
    end
end




"""

Parameter of a full order flux stator-rotor model without zero sequence flux in the stator.
The derivative of stator fluxes (ψd and ψq) is neglected. This is standard when
transmission network dynamics is neglected. Only one q-axis damping circuit
is considered. All per unit are in machine per unit.

Refer to Chapter 3 of Power System Stability and Control by P. Kundur or Chapter 11 of Power
System Dynamics: Stability and Control, by J. Machowski, J. Bialek and J. Bumby, for more details.
Note that the models are somewhat different (but equivalent) due to the different Park
Transformation used in both books.

# Conmutable structor
```julia
SimpleFullMachine(R, R_f, R_1d, R_1q, L_d, L_q, L_ad, L_aq, L_f1d, L_f, L_1d, L_1q)
```

#Arguments
* `R`::Float64 : Stator resistance after EMF in per unit
* `R_f`::Float64 : Field rotor winding resistance in per unit
* `R_1d`::Float64 : Damping rotor winding resistance on d-axis in per unit.
                  This value is denoted as RD in Machowski.
* `R_1q`::Float64 : Damping rotor winding resistance on q-axis in per unit.
                  This value is denoted as RQ in Machowski.
* `L_d`::Float64 : Inductance of fictitious damping that represent the effect
                  of the three-phase stator winding in the d-axis of the rotor, in per unit.
                  This value is denoted as L_ad + L_l in Kundur (and Ld in Machowski).
* `L_q`::Float64 : Inductance of fictitious damping that represent the effect
                  of the three-phase stator winding in the q-axis of the rotor, in per unit.
                  This value is denoted as L_aq + L_l in Kundur.
* `L_ad`::Float64 : Mutual inductance between stator winding and rotor field (and damping)
                    winding inductance on d-axis, in per unit
* `L_aq`::Float64 : Mutual inductance between stator winding and rotor damping winding
                    inductance on q-axis, in per unit
* `L_f1d`::Float64 : Mutual inductance between rotor field winding and rotor damping winding
                    inductance on d-axis, in per unit
* `L_ff`::Float64 : Field rotor winding inductance, in per unit
* `L_1d`::Float64 : Inductance of the d-axis rotor damping circuit, in per unit
* `L_1q`::Float64 : Inductance of the q-axis rotor damping circuit, in per unit
* `MVABase`::Float64 : Nominal Capacity in MVA

"""

mutable struct SimpleFullMachine <: Machine
  R::Float64
  R_f::Float64
  R_1d::Float64
  R_1q::Float64
  L_d::Float64
  L_q::Float64
  L_ad::Float64
  L_aq::Float64
  L_f1d::Float64
  L_ff::Float64
  L_1d::Float64
  L_1q::Float64
  MVABase::Float64
  inv_d_fluxlink::Array{Float64, 2}
  inv_q_fluxlink::Array{Float64, 2}
  n_states::Int64
  states::Vector{Symbol}
  ports::Ports

      function SimpleFullMachine(R::Float64,
                                R_f::Float64,
                                R_1d::Float64,
                                R_1q::Float64,
                                L_d::Float64,
                                L_q::Float64,
                                L_ad::Float64,
                                L_aq::Float64,
                                L_f1d::Float64,
                                L_ff::Float64,
                                L_1d::Float64,
                                L_1q::Float64,
                                MVABase::Float64)

            inv_d_fluxlink = inv( [ [-L_d L_ad L_ad]; #ψd: Eq 3.127 in Kundur
                                    [-L_ad L_ff L_f1d ]; #ψf: Eq 3.130 in Kundur
                                    [-L_ad L_f1d L_1d ] ]) #ψ1d: Eq 3.131 in Kundur
            inv_q_fluxlink = inv( [ [-L_q L_aq]; #ψq: Eq 3.128 in Kundur
                                    [-L_aq L_1q ] ]) #ψ1q: Eq 3.132 in Kundur
            n_states = 3
            states = [:ψf, :ψ1d, :ψ1q]
            @machine_ports

        new(R,
            R_f,
            R_1d,
            R_1q,
            L_d,
            L_q,
            L_ad,
            L_aq,
            L_f1d,
            L_ff,
            L_1d,
            L_1q,
            MVABase,
            inv_d_fluxlink,
            inv_q_fluxlink,
            n_states,
            states,
            Ports(state_input, inner_input))
    end
end
