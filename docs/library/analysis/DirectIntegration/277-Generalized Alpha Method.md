# Generalized Alpha

This command is used to construct a Generalized
$\alpha$ integration object. This is an implicit
method that like the HHT method allows for high frequency energy
dissipation and second order accuracy.
Depending on choices of input parameters, the method
can be unconditionally stable.

```tcl
integrator GeneralizedAlpha $alphaM $alphaF < $gamma $beta >
```

<hr />
<table>
<tbody>
<tr class="odd">
<td><p><code class="parameter-table-variable">alphaM</code></p></td>
<td><p>$\alpha_M$ factor</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">alphaF</code></p></td>
<td><p>$\alpha_F$ factor</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">gamma</code></p></td>
<td><p>$\gamma$ factor</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">beta</code></p></td>
<td><p>$\beta$ factor</p></td>
</tr>
</tbody>
</table>
<hr />

## Examples

```tcl
integrator GeneralizedAlpha 1.0 0.8
```
```python
{"integrator": ["GeneralizedAlpha", 1.0, 0.8]}
```

<p>NOTES:</p>
<ol>
<li>$\alpha_F$ and
$\alpha_M$ are defined differently that in the
paper, we use $\alpha_F = (1-\alpha_f)$ and
$\alpha_M=(1-\gamma_m)$ where
$\alpha_f$ and $\alpha_m$
are those used in the paper.</li>
<li>Like Newmark and all the implicit schemes, the unconditional
  stability of this method applies to linear problems. There are no
  results showing stability of this method over the wide range of
  nonlinear problems that potentially exist. Experience indicates that the
  time step for implicit schemes in nonlinear situations can be much
  greater than those for explicit schemes.</li>

<li>$\alpha_M = 1.0, \alpha_F = 1.0$ produces the Newmark Method.</li>
<li>$\alpha_M = 1.0$ corresponds to the HHT method.</li>
<li>The method is second-order accurate provided $\gamma = \tfrac{1}{2} + \alpha_M - \alpha_F$</li>

<li>The method is unconditionally stable provided $\alpha_M \ge \alpha_F \ge \tfrac{1}{2}, \beta \ge \tfrac{1}{4}
  +\tfrac{1}{2}(\gamma_M - \gamma_F)$
</li>
<li>$\gamma$ and $\beta$
  are optional. The default values ensure the method is unconditionally
  stable, second order accurate and high frequency dissipation is
  maximized.</li>
</ol>
<p>The defaults are:</p>
<dl>
<dt></dt>
<dd>

$$\gamma = \frac{1}{2} + \alpha_M - \alpha_F$$

</dd>
</dl>
<p>and</p>
<dl>
<dt></dt>
<dd>

$$\beta = \tfrac{1}{4}(1 + \gamma_M -
\gamma_F)^2$$

</dd>
</dl>

## Theory

<p>The Generalized $\alpha$ method (sometimes
called the $\alpha$ method) is a one step
implicit method for solving the transient problem which attempts to
increase the amount of numerical damping present without degrading the
order of accuracy. In the HHT method, the same Newmark approximations
are used:</p>
<dl>
<dt></dt>
<dd>

$$U_{t+\Delta t} = U_t + \Delta t \dot U_t + [(0.5 - \beta)
\Delta t^2] \ddot U_t + [\beta \Delta t^2] \ddot U_{t+\Delta
t}$$

</dd>
</dl>
<dl>
<dt></dt>
<dd>

$$\dot U_{t+\Delta t} = \dot U_t + [(1-\gamma)\Delta t] \ddot
U_t + [\gamma \Delta t ] \ddot U_{t+\Delta t} $$

</dd>
</dl>
<p>but the time-discrete momentum equation is modified:</p>

$$R_{t + \alpha_M \Delta t} = F_{t+\Delta t}^{ext} - M \ddot
U_{t + \alpha_M \Delta t} - C \dot U_{t+\alpha_F \Delta t} -
F^{int}(U_{t + \alpha_F \Delta t})
$$

<p>where the displacements and velocities at the intermediate point are
given by:</p>

$$U_{t+ \alpha_F \Delta t} = (1 - \alpha_F) U_t + \alpha_F
U_{t + \Delta t}$$

$$\dot U_{t+\alpha_F \Delta t} = (1-\alpha_F) \dot U_t +
\alpha_F \dot U_{t + \Delta t}$$

$$\ddot U_{t+\alpha_M \Delta t} = (1-\alpha_M) \ddot U_t +
\alpha_M \ddot U_{t + \Delta t}$$

<p>Following the methods outlined for Newmarks method, linearization of
the nonlinear momentum equation results in the following linear
equations:</p>

$$K_{t+\Delta t}^{*i} d U_{t+\Delta t}^{i+1} = R_{t+\Delta
t}^i$$

$$K_{t+\Delta t}^{*i} = \alpha_F K_t + \frac{\alpha_F
\gamma}{\beta \Delta t} C_t + \frac{\alpha_M}{\beta \Delta t^2}
M$$

<p>and</p>


$$R_{t+\Delta t}^i = F_{t + \Delta t}^{ext} - F(U_{t + \alpha
F \Delta t}^{i-1})^{int} - C \dot U_{t+\alpha F \Delta t}^{i-1} - M
\ddot U_{t+ \alpha M \Delta t}^{i-1}$$

<p>The linear equations are used to solve for 

$$U_{t+\alpha F
\Delta t}, \dot U_{t + \alpha F \Delta t} \ddot U_{t+ \alpha M \Delta
t}$$
Once convergence has been achieved the displacements,
velocities and accelerations at time $t + \Delta t$ can be computed.</p>

## References
<p>J. Chung, G.M.Hubert. "A Time Integration Algorithm for Structural
Dynamics with Improved Numerical Dissipation: The
Generalized-$\alpha$ Method" ASME Journal of
Applied Mechanics, 60, 371:375, 1993.</p>

<hr />

<p>Code Developed by: <span style="color:blue">fmk</span></p>

