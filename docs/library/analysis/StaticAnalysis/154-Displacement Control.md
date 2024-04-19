# Displacement Control

In an analysis step with Displacement Control we seek to
determine the time step $\Delta \lambda$ that will result in a displacement increment for
a particular degree-of-freedom at a node to be a prescribed value.

```tcl
integrator DisplacementControl $node $dof $incr
 < $numIter $DeltaUmin$ $DeltaUmax$ >
```

<hr />
<table>
<tbody>
<tr class="odd">
<td><p><code class="parameter-table-variable">node</code></p></td>
<td><p>node whose response controls solution</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">dof</code></p></td>
<td><p>degree of freedom at the node, valid options: 1 through <code>ndf</code> at node.</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">incr</code></p></td>
<td><p>first displacement increment $\Delta U_{\text{dof}}$ </p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">numIter</code></p></td>
<td><p>the number of iterations the user would like to occur in the
solution algorithm. Optional, default = 1.0.</p></td>
</tr>
<tr class="odd">
<td><p>$\Delta U\text{min}$</td>
<td><p>the min stepsize the user will allow. optional, defualt = $\Delta U_{min} = \Delta U_0$</p></td>
</tr>
<tr class="even">
<td><p>$\Delta U \text{max}$</p></td>
<td><p>the max stepsize the user will allow. optional, default = $\Delta U_{max} = \Delta U_0$</p></td>
</tr>
</tbody>
</table>

<hr />

## Examples

displacement control algorithm seking constant increment of 0.1 at node 1 at 2'nd dof.

```tcl
integrator DisplacementControl 1 2 0.1; 
```

```python
{"integrator": ["DisplacementControl", 1, 2, 0.1]}
```

<hr />

## Theory

If we write the governing finite element equation at $t + \Delta t$ as:

$$ 
R(\boldsymbol{u}_{n}, \lambda_{n}) = \lambda_{t+\Delta
t} \boldsymbol{p}_f - \boldsymbol{p}_{\sigma}(\boldsymbol{u}_{n})
$$


where \(\mathbf{f}_{\sigma}(\mathbf{u}_{n})\) are the internal
forces which are a function of the displacements \(\mathbf{u}_{n}\), \(\boldsymbol{p}_f\) is the set of
reference loads and \(\lambda\) is the load
multiplier. Linearizing the equation results in:


$$
\mathbf{K}_{n}^{*i} \Delta \mathbf{u}_{n}^{i+1} = \left(
\lambda^i_{n} + \Delta \lambda^i \right) \boldsymbol{p}_f -
\boldsymbol{p}_{\sigma}(\boldsymbol{u}_{n})
$$


This equation represents $n$ equations in $n+1$
unknowns, and so an additional equation is needed to solve the equation.
For displacement control, we introduce a new constraint equation in
which in each analysis step we set to ensure that the displacement
increment for the degree-of-freedom $\text{dof}$
at the specified node is:

$$
\Delta u_\text{dof} = \text{incr}
$$

### Incrementation

In Displacement Control $\Delta u_{\text{dof}}$ is incremented at $t +\Delta t$ to

$$
\Delta u_\text{dof}^{t+1} = \max \left( \Delta U_{\text{min}},
\min \left( \Delta U_{\text{max}},
\frac{\text{numIter}}{\text{lastNumIter}} \Delta U_\text{dof}^{t} \right) \right)
$$


<hr />

<p>Code Developed by: <span style="color:blue"> fmk
</span></p>
