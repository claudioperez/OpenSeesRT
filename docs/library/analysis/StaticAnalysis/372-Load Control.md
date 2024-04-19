# Load Control

This command is used to construct a LoadControl integrator
object.

```tcl
integrator LoadControl $lambda < $numIter $minLambda $maxLambda >
```
<hr />
<table>
<tbody>
<tr class="odd">
<td><p><code class="parameter-table-variable">lambda</code></p></td>
<td><p>the load factor increment $\lambda$ </p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">numIter</code></p></td>
<td><p>the number of iterations the user would like to occur in the
solution algorithm. Optional, default = 1.0.</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">minLambda</code></p></td>
<td><p>the min stepsize the user will allow. optional, defualt =
$\lambda_{\text{min}} = \lambda$</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">maxLambda</code></p></td>
<td><p>the max stepsize the user will allow. optional, default =
$\lambda_{\text{max}} = \lambda$</p></td>
</tr>
</tbody>
</table>
<hr />
<p>NOTES:</p>
<ol>
<li>The change in applied loads that this causes depends on the active
  load patterns (those load patterns not set constant) and the loads in
  the load patterns. If the only active loads acting on the domain are in
  load patterns with a Linear time series with a factor of 1.0, this
  integrator is the same as the classical load control method.</li>
<li>The optional arguments are supplied to speed up the step size in
  cases where convergence is too fast and slow down the step size in cases
  where convergence is too slow.</li>
</ol>
<hr />

## Examples

```tcl
integrator LoadControl 0.1;
```
```python
{"integrator": ["LoadControl", 0.1]};
```

<hr />

## Theory

### Incrementation
In Load Control the time in the domain is set to $\lambda + \lambda_{t+1}$ where,

$$
\Delta \lambda_i^1=\Delta \lambda_{i-1}^1\left(\frac{J_{\mathrm{d}}}{J_{i-1}}\right)
$$

$$\lambda_{t+1} = \max \left ( \lambda_{\text{min}}, \min \left (
\lambda_{\text{max}}, \frac{\text{numIter}}{\text{lastNumIter}} \lambda_{t}
\right ) \right ) $$

</dd>
</dl>
<hr />
<p>Code Developed by: <span style="color:blue">fmk</span></p>

# LoadControl 

```cpp
#include <analysis/integrator/LoadControl.h>
```


LoadControl is a subclass of StaticIntegrator, it is used to when
performing a static analysis on the `FE_Model` using the load control
method. In the load control method, the following constraint equation is
added to equation [staticFormTaylor](#staticFormTaylor) of the StaticIntegrator class:

$$
\lambda_n^{(i)} - \lambda_{n-1} = \delta \lambda_n
$$

where $\delta \lambda_n$ depends on $\delta \lambda_{n-1}$, the load
increment at the previous time step, $J_{n-1}$, the number of iterations
required to achieve convergence in the previous load step, and $J_d$, the
desired number of iteraions. $\delta
\lambda_n$ is bounded by $\delta \lambda_{\text{min}}$ and
$\delta \lambda_{\text{max}}$.
$$
\delta \lambda_n = \max \left( \delta \lambda_{\text{min}}, \min \left(
\frac{J_d}{J_{n-1}} \delta \lambda_{n-1}, \delta \lambda_{\text{max}} \right) \right)
$$

Knowing $\lambda_n^{(i)}$ prior to each iteration, the $N+1$ unknowns in
equation [staticFormTaylor](#staticFormTaylor), is reduced to $N$ unknowns and results in
the following equation:

$$
\boldsymbol{r}(\boldsymbol{u}_{n}) = \lambda_n^{(i)} \boldsymbol{p}_f 
 - \boldsymbol{p}_{\sigma}\left(\boldsymbol{u}_{n}^{(i)} \right) - 
\boldsymbol{K}_n^{(i)} 
(\boldsymbol{u}_{n} - \boldsymbol{u}_{n}^{(i)})  
\label{staticFormLoadControl}
$$

## Implementation

$\delta \lambda_1$ is the load factor used in the first step. The
arguments $J_d$, $\delta \lambda_{\text{min}}$, and $\delta
\lambda_{\text{max}}$ are used in the determination of the increment in the
load factor at each step.


```{.cpp}
int newStep(void);
```

The object obtains the current value of $\lambda$ from the AnalysisModel
object. It increments this by $\delta \lambda_n$.

$$
\delta \lambda_n = \max \left( \delta \lambda_{\text{min}}, \min \left(
\frac{J_d}{J_{n-1}} \delta \lambda_{n-1}, \delta \lambda_{\text{max}} \right) \right)
$$

It will then invoke `applyLoadDomain(0.0, $\lambda$)`{.cpp} on the
AnalysisModel object. Returns $0$ if successful. A warning message is
printed and a $-1$ is returned if no AnalysisModel is associated with
the object.

```{.cpp}
int update(const Vector &$\Delta \boldsymbol{u}$);
```

Invoked this causes the object to first increment the DOF_Group
displacements with $\Delta \boldsymbol{u}$, by invoking *incrDisp($\Delta \boldsymbol{u})$* on the
AnalysisModel, and then to update the domain, by invoking
`updateDomain()` on the AnalysisModel. Returns $0$ if successful, a
warning message and a $-1$ is returned if no AnalysisModel is associated
with the object.

Sets the value of the load increment in `newStep()` to be $\delta
\lambda$. Returns $0$.

```cpp
int sendSelf(int commitTag, Channel &theChannel);
```

Places in a vector if size 5 the value of $\delta \lambda_{n-1}$, $J_d$,
$J_{n-1}$, $\delta \lambda_{\text{min}}$ and $\delta \lambda_{\text{max}}$) and then
sends the Vector. Returns $0$ if successful, a warning message is
printed and a $-1$ is returned if `theChannel` fails to send the
Vector.

```cpp
int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
```

Receives in a Vector of size 5 the data that was sent in `sendSelf()`.
Returns $0$ if successful, a warning message is printed, $\delta \lambda$ 
is set to $0$, and a $-1$ is returned if *theChannel* fails to
receive the `Vector`.

```{.cpp}
int Print(OPS_Stream &s, int flag = 0);
```

The object sends to $s$ its type, the current value of $\lambda$, and
$\delta \lambda$.

