---
title: Static Analysis
description: |
  StaticAnalysis is a subclass of Analysis, it is used 
  to perform a static analysis on the Domain.
...

# Static Analysis

In static nonlinear finite element problems we seek a solution ($\boldsymbol{u}$, $\lambda$) 
to the nonlinear vector function

$$
\boldsymbol{r}(\boldsymbol{u}, \lambda) = \lambda \boldsymbol{p}_f - \boldsymbol{p}_r(\boldsymbol{u}) = \boldsymbol{0}
%\label{staticGenForm}
$$

The most widely used technique for solving the non-linear finite element
equation, equation [femGenForm](#femGenForm), is to use 
an incremental scheme. In the incremental formulation, a solution to 
the equation is sought at successive incremental steps.

$$
\boldsymbol{r}(\boldsymbol{u}_{n}, \lambda_n) = \lambda_n \boldsymbol{p}_f - \boldsymbol{p}_{\sigma}(\boldsymbol{u}_{n})
$$

The solution of this equation is typically obtained using an iterative
procedure, in which a sequence of approximations ($\boldsymbol{u}_{n}^{(i)}$,
$\lambda_n^{(i)}$), $i=1,2, ..$ is obtained which converges to the
solution ($\boldsymbol{u}_n$, $\lambda_n)$. The most frequently used iterative
schemes, such as Newton-Raphson, modified Newton, and quasi Newton
schemes, are based on a Taylor expansion of
equation [staticIncForm](#staticIncForm) about ($\boldsymbol{u}_{n}$, $\lambda_n$):

$$
\boldsymbol{r}(\boldsymbol{u}, \lambda)(\boldsymbol{u}_{n},\lambda_n) = \lambda_n^{(i)} \boldsymbol{p}_f 
 - \boldsymbol{p}_{\sigma}\left(\boldsymbol{u}_{n}^{(i)} \right) - \left[
\begin{array}{cc}
\boldsymbol{K}_n^{(i)} & -\boldsymbol{p}_f \\
\end{array} \right] 
\left\{
\begin{array}{c}
\boldsymbol{u}_{n} - \boldsymbol{u}_{n}^{(i)}  \\ 
\lambda_n - \lambda_n^{(i)} 
\end{array} \right\}
%\label{staticFormTaylor}
$$

which a system of of $N$ equations with ($N+1$) unknowns. Two solve
this, an additional equation is required, the constraint equation. The
constraint equation used depends on the static integration scheme, of
which there are a number, for example load control, arc length, and
displacement control.


The following are the aggregates of such an analysis type:

-  **StaticIntegrator** - an algorithmic class which provides methods
   which are invoked by the `FE_Element` to determine their current
   tangent and residual matrices; that is this is the class that sets
   up the system of equations. It also provides the `update()` method
   which is invoked to set up the appropriate dof response values once
   the solution algorithm has formed and solved the system of
   equations.
   - [Load Control]() -- Specifies the incremental load factor to be applied to the loads in the domain
   - [Displacement Control]() -- Specifies the incremental displacement at a specified DOF in the domain
   - [Minimum Unbalanced Displacement Norm](StaticIntegrator/) -- Specifies the incremental load factor such that the residual displacement norm in minimized
   - [Arc Length](StaticIntegrator/ArcLength1) -- Specifies the incremental arc-length of the load-displacement path


-  **EquiSolnAlgo** - an algorithmic class specifying the sequence of
   operations to be performed in setting up and solving the finite
   element equation which can be represented by the equation 
   $$
   \boldsymbol{K}_{\mathrm{eff}}(\boldsymbol{u}) \boldsymbol{u} = \boldsymbol{p}(\boldsymbol{u})
   $$

-  **ConstraintHandler** - a class which creates the `DOF_Group` and
   `FE_Element` objects, the type of objects created depending on how the
   specified constraints in the domain are to be handled.

-  **DOF_Numberer** - a class responsible for providing equation
   numbers to the individual degrees of freedom in each DOF_Group
   object.

-  **LinearSOE** - a numeric class responsible for the creation and
   subsequent solution of large systems of linear equations of the form
   $Ax = b$, where $A$ is a matrix and $x$ and $b$ are vectors.



|                | `-det`   |  `-iterScale` |  `-exp`  |  `phat`  | `-min/max` |
|----------------|----------|---------------|----------|----------|------|
| `MinUnbalDisp` | Fix      | 
| `ArcLength`    | Add      | Add           |   Add    |
| `LoadControl`  | Add      | 
| `DispControl`  | Add      | 



|              | `-det`  |
|--------------|---------|
| `Umfpack`    | 
| `FullGen`    |
| `BandGen`    |
| `ProfileSPD` |
| `Band SPD`   |
| `SparseGen`  |
| `SparseSym`  |
| `UmfGen`     |
| `Conjugate Gradient`
