" Define a step change of the mechanical power at t = 1.0"

function change_t_one(x,t,integrator)
  t in [1.0]
end

function Y_change!(integrator)
  integrator.p[2].Ybus = integrator.p[1]
end
