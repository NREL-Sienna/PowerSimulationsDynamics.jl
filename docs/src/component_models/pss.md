# Power System Stabilizers (PSS)

PSS are used to add an additional signal ``v_s`` to the field voltage: ``v_f = v_f^{\text{avr}} + v_s``.

## Fixed PSS ```[PSSFixed]```

This is a simple model that set the stabilization signal to be equal to a desired constant value ``v_s = v_{s}^{\text{fix}}``. The absence of PSS can be modelled using this component with ``v_s^{\text{fix}} = 0``.

## Simple PSS ```[PSSSimple]```

This is the most basic PSS that can be implemented, on which the stabilization signal is  a proportional controller over the frequency and electrical power:

```math
\begin{align}
v_s = K_{\omega}(\omega - \omega_s) + K_p(\omega \tau_e - P_{\text{ref}}) \tag{12a}
\end{align}
```
