# Central Difference

<p>This command is used to construct a Central Difference integrator
object.</p>

```tcl
integrator CentralDifference
```
<hr />

## Examples

```tcl
integrator CentralDifference
```

```python
{"integrator": ["CentralDifference"]}
```

<hr />

<p>NOTES:</p>

- The calculation of $U_{t+\Delta t}$, as shown
  below, is based on using the equilibrium equation at time t. For this
  reason the method is called an <strong>explicit integration
  method</strong>.
- If there is no rayleigh damping and the $C$ matrix is 0, for a
  diagonal mass matrix a diagonal solver may and should be used.
- For stability, $\frac{\Delta t}{T_n} < \frac{1}{\pi}$

<hr />

## Theory

The Central difference approximations for velocity and
acceleration:

$$v_n = \frac{u_{n+1} - u_{n-1}}{2 \Delta t}$$



$$a_n = \frac{u_{n+1} - 2 u_n + u_{n-1}}{\Delta
t^2}$$


In the Central Difference method we determine the displacement
solution at time $t+\Delta t$ by considering the
the eqilibrium equation for the finite element system in motion at time
$t$:

$$M \ddot{\boldsymbol{u}}_t + C \dot{\boldsymbol{u}}_t + K {\boldsymbol{u}}_t = R_t$$

which when using the above two expressions of becomes:

$$\left ( \frac{1}{\Delta t^2} M + \frac{1}{2 \Delta t} C
\right ) U_{t+\Delta t} = R_t - \left (K - \frac{2}{\Delta t^2}M \right
)U_t - \left (\frac{1}{\Delta t^2}M - \frac{1}{2 \Delta t} C \right)
U_{t-\Delta t} $$



## References

<hr />

Code Developed by: <span style="color:blue"> fmk
</span>
