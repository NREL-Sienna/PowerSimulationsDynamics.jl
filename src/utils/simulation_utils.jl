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
  integrator.p[2].dyn_injections[1].P_ref += integrator.p[1][1]
end

function step_change4!(integrator)
  integrator.p[2].dyn_injections[2].P_ref += integrator.p[1][1]
end

function Y_change!(integrator)
  integrator.p[2].Ybus = integrator.p[1]
end
