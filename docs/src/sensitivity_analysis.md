# Sensitivity Analysis 

PowerSimulationsDynamics has limited support for performing sensitivity analysis of system parameters with respect to a user defined loss function. See the tutorial for an example of sensitivity analysis. 

* `ForwardDiffSensitivity()` is used as the method for differentiating through the solve. See `SciMLSensitivity.jl` for more details. 
* The gradient of the function provided by PSID can be calculated using `Zygote.jl`.
* The Jacobian is not passed to the ODE Problem during sensitivity analysis. 
* Parameters for sensitivity analysis must not change the steady state operating condition. Parameters cannot be in the mass matrix. This limitation is expected to be relaxed in the future. 
* Limited to mass matrix formulations and pure Julia solvers. Can check if a solver is compatible with `SciMLBase.isautodifferentiable()`
