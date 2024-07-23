# Minimum Unbalanced Displacement Norm

```tcl
integrator MinUnbalDispNorm $dlambda1 < $Jd $minLambda $maxLambda >
```

<hr />
<table>
<tbody>
<tr class="odd">
<td><p><code class="parameter-table-variable">dlambda1</code></p></td>
<td><p>first load increment (pseudo-time step) at the first iteration in
the next invocation of the <a href="Analysis_Command"
title="wikilink">analysis</a> command.</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">Jd</code></p></td>
<td><p>factor relating first load increment at subsequent time steps
(optional, default: 1.0)</p></td>
</tr>
<tr class="odd">
<td><p><strong>$minLambda, $maxLambda</strong></p></td>
<td><p>arguments used to bound the load increment (optional, default:
$dLambda1 for both)</p></td>
</tr>
</tbody>
</table>
<hr />

## Examples

```tcl
integrator MinUnbalDispNorm 0.1;
```

```python
{"integrator": ["MinUnbalDispNorm", 0.1]};
```

<hr />

## Theory

### Continuation

In this instance, the constraint equation involving $\Delta \lambda_i^j$ is
$$
\frac{\partial}{\partial \Delta \lambda_i^j}\left[\{\Delta u\}_i^{j} \cdot \{\Delta \delta\}_i^j\right]=0
$$

which guarantees a minimum value for the unbalanced displacement norm in each iteration. Expanding $\{\Delta \delta\}_i^j$ as defined in equation (9) and evaluating equation (28) furnishes
$$
\Delta \lambda_i^j=\frac{-\left\{\delta_1\right\}_i^{\mathrm{T}}\left\{\Delta \delta_{\mathrm{R}}\right\}_i^j}{\left\{\delta_1\right\}_i^{\mathrm{T}}\left\{\delta_{\mathrm{I}}\right\}_i}
$$

### Incrementation

The load increment at iteration $i$, $d\lambda_{1,i}$ , is related to the load
increment at $(i-1)$, $d\lambda_{1,i-1}$, and the number of iterations at
$(i-1)$, $J_{i-1}$, by the following:


$$
d\lambda_{1,i} = d\lambda_{1,i-1}\frac{J_d}{J_{i-1}}
$$

<hr />
<p>Code Developed by: <span style="color:blue"> fmk
</span></p>
