function step_change_first!(integrator)
  integrator.p[2].dyn_injections[1].P_ref += integrator.p[1][1]
end

function step_change_second!(integrator)
    integrator.p[2].dyn_injections[2].P_ref += integrator.p[1][1]
end
