" Define a step change of the mechanical power at t = 1.0"


function change_t_one(x,t,integrator)
  t in [1.0]
end

function step_change!(integrator)
  integrator.u0[1] = integrator.u0[1] + 0.3
end

function step_change1!(integrator)
  integrator.p[2].dyn_injections[1].tg.P_order += 0.1
end

function step_change2!(integrator)
  integrator.OMIB.dyn_injections[1].tg.P_order = integrator.OMIB.dyn_injections[1].tg.P_order + 0.3
end

function step_change3!(integrator)
  integrator.p[2].dyn_injections[1].P_ref += 0.1
end

function Y_change!(integrator)
  integrator.p[2].Ybus = Ybus_fault
end


"""

Function to obtain which indices of the vector of the device states correspond
to the states of the component (of that device) selected.

The function receives a dictionary with the symbol and position (dev_ix) of
the device (e.g. the generator) and also receives the specific component of
that device (e.g. AVR, machine, PSS, Governor).

"""


function _get_index(dev_ix::Dict{Symbol, Int64},
                   gen_comp::G) where {G <: GeneratorComponent}

  states = gen_comp.states
  n_state =
  comp_ix = Vector{Int64}(undef, length(states))
  for (ix, s) in enumerate(states)
    comp_ix[ix] = get(dev_ix, s, -99)
    #look into base.pairs for performance improvement
  end
  return comp_ix
end

function get_states(component::C,
                    states,
                    device_ix::Dict{Symbol,Int64},) where {C <: GeneratorComponent}



end
"""

Function to obtain AVR states.
The function only returns the fixed EMF in case of AVRFixed type.
In the other cases returns the vector of states of the AVR.

"""

function get_avr_states(avr::A, states, gen_ix::Dict{Symbol,Int64}) where  {A <: AVR}
  avr_x = typeof(avr) == AVRFixed ? avr.Emf : view(states, get_index(gen_ix, avr))
  return avr_x
end

"""

Function to obtain PSS states.
The function only returns the fixed V_ref in case of PSSFixed type.
In the other cases returns the vector of states of the PSS.

"""

function get_pss_states(pss::P, states, gen_ix::Dict{Symbol,Int64}) where {P <: PSS}
    pss_x = typeof(pss) == PSSFixed ? pss.V_ref : view(states, get_index(gen_ix, pss))
  return pss_x
end

"""

Function to obtain Turbine Governor states.
The function returns an empty vector in case of TGFixed type.
In the other cases returns the vector of states of the TG.

"""

function get_tg_states(tg::TG, states, gen_ix::Dict{Symbol,Int64}) where {TG <: TurbineGov}
    tg_x = typeof(tg) == TGFixed ? Vector{Float64}() : view(states, get_index(gen_ix, tg))
  return tg_x
end

#function get_index(device_index::Dict{Symbol, Int64}, component::C) where {C <: GeneratorComponent}
#    return sort([v for (k,v) in device_index if k in component.states])
#end
