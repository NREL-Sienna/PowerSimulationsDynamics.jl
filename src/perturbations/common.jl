#=
function change_t_one(x,t,integrator)
  t in [1.0]
end

function Y_change_ext!(integrator)
  PSY.get_ext(integrator.p[2])["Ybus"] = integrator.p[1]
end
cb = DiffEqBase.DiscreteCallback(change_t_one, Y_change_ext!)
=#
