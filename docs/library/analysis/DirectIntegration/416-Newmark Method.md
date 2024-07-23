# Newmark

This command is used to construct a Newmark integrator object.

:::{apidoc="opensees.library.analysis.Newmark"}

```tcl
integrator Newmark $gamma $beta
```
<hr />
<table>
<tbody>
<tr class="odd">
<td><code class="parameter-table-variable">gamma</code></td>
<td>$\gamma$ factor</td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">beta</code></td>
<td>$\beta$ factor</td>
</tr>
</tbody>
</table>
<hr />

:::

## Examples

```tcl
integrator Newmark 0.5 0.25
```

```python
model.integrator("Newmark", beta=0.25, gamma=0.5, solve="a")
```

<hr />
NOTES:

- If the accelerations are chosen as the unknowns and
  $\beta$ is chosen as $0$, the formulation results
  in the fast but conditionally stable explicit Central Difference method.
  Otherwise the method is implicit and requires an iterative solution
  process.

- Two common sets of choices are

  - Average Acceleration Method ($\gamma=\frac{1}{2}, \beta = \frac{1}{4}$
  - Linear Acceleration Method ($\gamma=\frac{1}{2}, \beta = \frac{1}{6}$

- $\gamma > \frac{1}{2}$ results in numerical damping proportional to $\gamma - \frac{1}{2}$
- The method is second order accurate if and only if $\gamma=\frac{1}{2}$
- The method is conditionally stable for $\beta \ge \frac{\gamma}{2} \ge \frac{1}{4}$

<hr />

## Theory

The Newmark method is a one step implicit method for solving the
transient problem, represented by the residual for the momentum
equation:
$$
R_{t + \Delta t} = F_{t+\Delta t}^{ext} - M \ddot{\boldsymbol{u}}_{t +
\Delta t} - C \dot{\boldsymbol{u}}_{t + \Delta t} + F(U_{t + \Delta
t})^{int}
$$



Using the Taylor series approximation of $U_{t+\Delta t}$ and $\dot{\boldsymbol{u}}_{t+\Delta t}$:
$$
\boldsymbol{u}_{t+\Delta t} = \boldsymbol{u}_t + \Delta t \dot{\boldsymbol{u}}_t + \frac{\Delta
t^2}{2} \ddot{\boldsymbol{u}}_t + \frac{\Delta t^3}{6} \ddot{\boldsymbol{u}}_t + \cdots
$$

$$\dot{\boldsymbol{u}}_{t+\Delta t} = \dot{\boldsymbol{u}}_t + \Delta t \ddot{\boldsymbol{u}}_t +
\frac{\Delta t^2}{2} \ddot{\boldsymbol{u}}_t + \cdots 
$$

Newton truncated these using the following:
$$
\boldsymbol{u}_{t+\Delta t} = \boldsymbol{u}_t + \Delta t \dot{\boldsymbol{u}}_t + \frac{\Delta t^2}{2} \ddot{\boldsymbol{u}} + \beta {\Delta t^3} \ddot{\boldsymbol{u}}
$$

$$
\dot{\boldsymbol{u}}_{t + \Delta t} = \dot{\boldsymbol{u}}_t + \Delta t \ddot{\boldsymbol{u}}_t +
\gamma \Delta t^2  \ddot{\boldsymbol{u}}
$$

in which he assumed linear acceleration within a time step, i.e.

$$
\ddot{\boldsymbol{u}} = \frac{{\ddot{\boldsymbol{u}}_{t+\Delta t}} - \ddot{\boldsymbol{u}}_t}{\Delta t}
$$
which results in the following expressions:

$$
U_{t+\Delta t} = U_t + \Delta t \dot{\boldsymbol{u}}_t + [(0.5 - \beta)
\Delta t^2] \ddot{\boldsymbol{u}}_t + [\beta \Delta t^2] \ddot{\boldsymbol{u}}_{t+\Delta t}
$$

$$\dot{\boldsymbol{u}}_{t+\Delta t} = \dot{\boldsymbol{u}}_t + [(1-\gamma)\Delta t] \ddot{\boldsymbol{u}}_t + [\gamma \Delta t ] \ddot{\boldsymbol{u}}_{t+\Delta t} 
$$

The variables $\beta$ and
$\gamma$ are numerical parameters that control
both the stability of the method and the amount of numerical damping
introduced into the system by the method. For
$\gamma=\frac{1}{2}$ there is no numerical
damping; for $\gamma \ge \frac{1}{2}$ numerical
damping is introduced.

The linearization of the Newmark equations gives:


$$dU_{t+\Delta t}^{i+1} = \beta \Delta t^2 d \ddot{\boldsymbol{u}}_{t+\Delta t}^{i+1}$$



$$d \dot{\boldsymbol{u}}_{t+\Delta t}^{i+1} = \gamma \Delta t \ddot{\boldsymbol{u}}_{t+\Delta t}^{i+1}$$


which gives the update formula when displacement increment is used as
unknown in the linearized system as:


$$U_{t+\Delta t}^{i+1} = U_{t+\Delta t}^i + dU_{t+\Delta t}^{i+1}$$



$$
\dot{\boldsymbol{u}}_{t+\Delta t}^{i+1} 
= \dot{\boldsymbol{u}}_{t+\Delta t}^i 
+ \frac{\gamma}{\beta \Delta t} dU_{t+\Delta t}^{i+1}
$$



$$
\ddot{\boldsymbol{u}}_{t+\Delta t}^{i+1} 
= \ddot{\boldsymbol{u}}_{t+\Delta t}^i 
+ \frac{1}{\beta \Delta t^2}dU_{t+\Delta t}^{i+1}
$$


The linearization of the momentum equation using the displacements as
the unknowns leads to the following linear equation:


$$
K_{t+\Delta t}^{*i} \Delta U_{t+\Delta t}^{i+1} = R_{t+\Delta t}^i
$$

where


$$
K_{t+\Delta t}^{*i} = K_t + \frac{\gamma}{\beta \Delta t}
C_t + \frac{1}{\beta \Delta t^2} M
$$

and


$$
R_{t+\Delta t}^i = F_{t + \Delta t}^{ext} - F(U_{t + \Delta
t}^{i-1})^{int} - C \dot{\boldsymbol{u}}_{t+\Delta t}^{i-1} - M \ddot{\boldsymbol{u}}_{t+ \Delta
t}^{i-1}
$$


## References
Newmark, N.M. "A Method of Computation for Structural Dynamics" ASCE
Journal of Engineering Mechanics Division, Vol 85. No EM3, 1959.

<hr />
Code Developed by: <span style="color:blue"> fmk
</span>
