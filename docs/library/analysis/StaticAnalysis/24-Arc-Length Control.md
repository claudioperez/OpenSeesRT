# Arc-Length Control

This command is used to construct an ArcLength integrator object. In
an analysis step with ArcLength we seek to determine the time step that
will result in our constraint equation being satisfied.

```tcl
integrator ArcLength $s $alpha
```

<hr />

<table>
<tbody>
<tr class="odd">
<td><p><code class="parameter-table-variable">s</code></p></td>
<td><p>$s$ the arcLength.</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">alpha</code></p></td>
<td><p>$\alpha$ a scaling factor on the reference
loads.</p></td>
</tr>
</tbody>
</table>

<hr />

## Examples

```{.tcl .example}
integrator ArcLength 1.0 0.1;
```

```{.py .example}
{"integrator": ["ArcLength", 1.0, 0.1]}
```

<hr />

## Theory

If we write the governing finite element equation at $t + \Delta t$ as:

$$ 
\boldsymbol{r}(\boldsymbol{u}_{n}, \lambda_{n}) 
  = \lambda_{n} \boldsymbol{p}_f - \boldsymbol{p}_{\sigma}(\boldsymbol{u}_{n})
$$

where $\boldsymbol{p}_{\sigma}(\boldsymbol{u}_{n})$ are the internal
forces which are a function of the displacements $\boldsymbol{u}_{n}$,
$\boldsymbol{p}_f$ is the set of reference applied loads and $\lambda$ is the load
multiplier. Linearizing the equation results in:

$$
\boldsymbol{K}_{n}^{*i} \Delta \boldsymbol{u}_{n}^{i+1} = \left(
\lambda^i_{n} + \Delta \lambda^i \right) \boldsymbol{p}_f -
\boldsymbol{p}_{\sigma}(\boldsymbol{u}_{n})
$$

This equation represents $n$ equations in $n+1$
unknowns, and so an additional equation is needed to solve the equation.
For displacement control, we introduce a new constraint equation in
which in each analysis step we set to ensure that the displacement
increment for the degree-of-freedom `dof`
at the specified node is:

<hr />

On the first equilibrium iteration within a load step, both arc length integrators compute the change in load factor, $\Delta \lambda^1$, from the "tangent displacements", $\mathbf{U}_{f T}=\mathbf{K}_T \backslash \mathbf{P}_{r e f}$ , and the user specified arc length and $\alpha$.
$$
\Delta \lambda^1=\sqrt{\frac{s^2}{\mathbf{U}_{f T} \cdot \mathbf{U}_{f T}+\alpha^2}}
$$

<p>Code Developed by: <span style="color:blue"> fmk</span></p>

