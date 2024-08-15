# DirectIntegrationAnalysis 

```python
class DirectIntegrationAnalysis(model, patterns, strategy)
```

`DirectIntegrationAnalysis` is used to perform a transient analysis using an
incremental approach on the `model`.


The following are the aggregates of such an
analysis type:

-   `model` (**AnalysisModel**) - a container class holding the `FE_Element` and
    `DOF_Group` objects created by the ConstraintHandler object.

-   `integrator` (**TransientIntegrator**) - an algorithmic class which provides
    methods which are invoked by the `FE_Element` to determine their
    current tangent and residual matrices. That is, this is the class
    that sets up the system of equations. 

    $$\mathbf{M}\ddot{\boldsymbol{u}}+\mathbf{C}\dot{\boldsymbol{u}}+\boldsymbol{p}_{\sigma}(\boldsymbol{u})=\boldsymbol{p}_f(t)$$

    It also provides the `commit()` method which is invoked to set up the
    appropriate dof response values once the solution algorithm has formed and
    solved the system of equations.

-   `constraints` (**ConstraintHandler**) - a class which creates the `DOF_Group` and
    `FE_Element` objects, the type of objects created depending on how the
    specified constraints in the domain are to be handled.

-   `numberer` (**DOF_Numberer**) - a class responsible for providing equation
    numbers to the individual degrees of freedom in each `DOF_Group`
    object.

-   `system` (**LinearSOE**) - a numeric class responsible for the creation and
    subsequent solution of large systems of linear equations of the form
    $Ax = b$, where $A$ is a matrix and $x$ and $b$ are vectors.

-   `algorithm` (**EquiSolnAlgo**) - an algorithmic class specifying the sequence of
    operations to be performed in setting up and solving the finite
    element equation which can be represented by the equation $K(\boldsymbol{u}) \boldsymbol{u} = P(\boldsymbol{u})$.


### Transient Integrators
Determing the next time step for an analysis including inertial effects is done by the following schemes

- [Newmark](TransientIntegrator/Newmark) -- The two parameter time-stepping method developed by Newmark
- [HHT](TransientIntegrator/HHT) -- The three parameter Hilbert-Hughes-Taylor time-stepping method
- [Generalized Alpha](TransientIntegrator/GeneralizedAlpha) -- Generalization of the HHT algorithm with improved numerical damping
- [Central Difference]() -- Approximates velocity and acceleration by centered finite differences of displacement
- [TRBDF2](TRBDF2) --  A composite scheme that alternates between the Trapezoidal scheme and a 3 point backward Euler scheme.
- [Explicit difference](Explicitdifference) -- 

## Theory

In nonlinear transient finite element problems we seek a solution 
($\boldsymbol{u}$, $\dot{\boldsymbol{u}}$, $\ddot{\boldsymbol{u}}$) to the nonlinear vector equation

$$\boldsymbol{R}({\boldsymbol{u}},\dot{\boldsymbol{u}}, \ddot{\boldsymbol{u}}) = \boldsymbol{p}_f(t) - \boldsymbol{p}_{\mathrm{i}}(\ddot{\boldsymbol{u}}) - \boldsymbol{p}_{\sigma}({\boldsymbol{u}}, \dot{\boldsymbol{u}}) = \boldsymbol{0}$$
{#femGenForm}

The most widely used technique for solving the transient non-linear finite element equation,
equation [\[femGenForm](#femGenForm){reference-type="ref"
reference="femGenForm"}, is to use an incremental direct integration
scheme. In the incremental formulation, a solution to the equation is
sought at successive time steps $\Delta t$ apart.

$$\boldsymbol{R}({\boldsymbol{u}}_{n \Delta t},\dot{\boldsymbol{u}}_{n \Delta t}, \ddot{\boldsymbol{u}}_{n \Delta t}) = \boldsymbol{p}_f(n \Delta t) -
\boldsymbol{p}_{\mathrm{i}}(\ddot{\boldsymbol{u}}_{n \Delta t}) - \boldsymbol{p}_{\sigma}({\boldsymbol{u}}_{n \Delta t}, \dot{\boldsymbol{u}}_{n \Delta t})
$$
{#fullTimeForm}

For each time step, $t$, the integration schemes provide two operators,
$\operatorname{I}_1$ and $\operatorname{I}_2$, to relate the velocity and accelerations at the time
step as a function of the displacement at the time step and the response
at previous time steps:

$$\dot {\boldsymbol{u}}_{t} = {\mathrm{I}}_1 ({\boldsymbol{u}}_t, {\boldsymbol{u}}_{t-\Delta t}, \dot {\boldsymbol{u}}_{t-\Delta t},
\ddot {\boldsymbol{u}}_{t - \Delta t}, {\boldsymbol{u}}_{t - 2\Delta t}, \dot {\boldsymbol{u}}_{t - 2 \Delta t}. ..., )
%\label{I1}
$$

$$\ddot {\boldsymbol{u}}_{t} = {\mathrm{I}}_2 ({\boldsymbol{u}}_t, {\boldsymbol{u}}_{t-\Delta t}, \dot{\boldsymbol{u}}_{t-\Delta t},
\ddot{\boldsymbol{u}}_{t - \Delta t}, {\boldsymbol{u}}_{t - 2\Delta t}, \dot{\boldsymbol{u}}_{t - 2 \Delta t}. ..., )
%\label{I2}
$$

These allow us to rewrite
equation [\[fullTimeForm](#fullTimeForm){reference-type="ref"
reference="fullTimeForm"}, in terms of a single response quantity,
typically the displacement:

$$
\boldsymbol{r}({\boldsymbol{u}}_t) = \boldsymbol{p}_f(t) - \boldsymbol{p}_{\mathrm{i}}(\ddot{\boldsymbol{u}}_t) - \boldsymbol{p}_{\sigma}({\boldsymbol{u}}_t, \dot{\boldsymbol{u}}_t)
%\label{genForm}
$$

The solution of this equation is typically obtained using an iterative
procedure, i.e. making an initial prediction for ${\boldsymbol{u}}_{t}$, denoted
${\boldsymbol{u}}_{t}^{(0)}$ a sequence of approximations ${\boldsymbol{u}}_{t}^{(i)}$, $i=1,2, ..$
is obtained which converges (we hope) to the solution ${\boldsymbol{u}}_{t}$. The most
frequently used iterative schemes, such as Newton-Raphson, modified
Newton, and quasi Newton schemes, are based on a Taylor expansion of
equation [genForm](#genForm) about ${\boldsymbol{u}}_{t}$:

$$\boldsymbol{r}({\boldsymbol{u}}_{t}) = 
\boldsymbol{r}({\boldsymbol{u}}_{t}^{(i)}) +
\left[ {\frac{\partial \boldsymbol{r}}{\partial {\boldsymbol{u}}_t} \vert}_{{\boldsymbol{u}}_{t}^{(i)}}\right]
\left( {\boldsymbol{u}}_{t} - {\boldsymbol{u}}_{t}^{(i)} \right)$$

$$
\boldsymbol{r}({\boldsymbol{u}}_{t}) = \boldsymbol{p}_f (t) - \boldsymbol{p}_{\mathrm{i}} \left( \ddot {\boldsymbol{u}}_{t}^{(i)} \right) - \boldsymbol{p}_{\sigma} \left( \dot {\boldsymbol{u}}_{t}^{(i)}, {\boldsymbol{u}}_{t}^{(i)} \right)- \left[
  \boldsymbol{M}^{(i)} {\mathrm{I}}_2'
+ \boldsymbol{C}^{(i)} {\mathrm{I}}_1'
+ \boldsymbol{K}^{(i)}  \right]
 \left( {\boldsymbol{u}}_{t} - {\boldsymbol{u}}_{t}^{(i)} \right)
%\label{femGenFormTaylor}
$$

To start the iteration scheme, trial values for ${\boldsymbol{u}}_{t}$, $\dot
{\boldsymbol{u}}_{t}$ and $\ddot {\boldsymbol{u}}_{t}$ are required. These are obtained by assuming
${\boldsymbol{u}}_{t}^{(0)} = {\boldsymbol{u}}_{t-\Delta t}$. The $\dot {\boldsymbol{u}}_{t}^{(0)}$ and
$\ddot {\boldsymbol{u}}_{t}^{(0)}$ can then be obtained from the operators for the
integration scheme.
Subclasses of `TransientIntegrator` provide methods informing the
`FE_Element` and `DOF_Group` objects how to build the tangent and residual
matrices and vectors. They also provide the method for updating the
response quantities at the DOFs with appropriate values; these values
being some function of the solution to the linear system of equations.



## Implementation

```cpp
#include <analysis/analysis/DirectIntegrationAnalysis.h>

class DirectIntegrationAnalysis: public TransientAnalysis;
```

### Constructor and Destructor

```cpp
DirectIntegrationAnalysis::DirectIntegrationAnalysis(
                 Domain              &the_Domain,
                 ConstraintHandler   &theHandler,
                 DOF_Numberer        &theNumberer,
                 AnalysisModel       &theModel,
                 EquiSolnAlgo        &theSolnAlgo,
                 LinearSOE           &theLinSOE,
                 TransientIntegrator &theTransientIntegrator,
                 ConvergenceTest     *theConvergenceTest,
                 int                 num_SubLevels,
                 int                 num_SubSteps
)
```


The constructor is responsible for setting up all
links needed by the objects in the aggregation. It invokes
`setLinks(theDomain)` on `theModel`,
`setLinks(theDomain,theModel,theIntegrator)` on `theHandler`,
`setLinks(theModel)` on `theNumberer`, `setLinks(theModel, theSOE)`{.cpp} on
`theIntegrator` and `setLinks(theModel,theAnalysis, theIntegrator, theSOE)` 
on `theSolnAlgo`.
Sets theModel and theSysOFEqn to 0 and the Algorithm to the one supplied.

```cpp
~DirectIntegrationAnalysis
```
Does nothing. `clearAll()` must be invoked if the destructor on the
objects in the aggregation need to be invoked.

### Public Methods

:::{.admonition}
```cpp
int analyze(int steps, double dT);
```
:::
Invokes `analyzeStep(dT)` `steps` number of times.

:::{.admonition}
```cpp
int analyzeStep(double dT);
```
:::
Invoked to perform a transient analysis on the `FE_Model`. The method
checks to see if the domain has changed before it performs the analysis.
The `DirectIntegrationAnalysis` object performs the following:

```cpp
theAnalysisModel->analysisStep(dT);

if (theDomain->hasDomainChanged() != this->domainStamp)
  this->domainChanged();

theIntegrator->newStep(dT);
theAlgorithm->solveCurrentStep();
theIntegrator->commit();
```

The type of analysis performed, depends on the type of the objects in
the analysis aggregation. If any of the methods invoked returns a
negative number, an error message is printed, `revertToLastCommit()` is
invoked on the Domain, and a negative number is immediately returned.
Returns a $0$ if the algorithm is successful.

:::{.admonition}
```{.cpp}
void clearAll(void);
```
:::
Will invoke the destructor on all the objects in the aggregation. NOTE
this means they must have been constructed using `new()`, otherwise a
segmentation fault can occur.

:::{.admonition}
```{.cpp}
void domainChange(void);
```
:::
This is a method invoked by a domain which indicates to the analysis
that the domain has changed. The method invokes the following:

1. It invokes `clearAll()` on `theModel` which causes the `AnalysisModel`
   to clear out its list of `FE_Elements` and `DOF_Groups`, and
   `clearAll()` on `theHandler`.

2. It then invokes `handle()` on `theHandler`. This causes the
   constraint handler to recreate the appropriate `FE_Element` and
   `DOF_Groups` to perform the analysis subject to the boundary
   conditions in the modified domain.

3. It then invokes `numberDOF()` on `theNumberer`. This causes the DOF
   numberer to assign equation numbers to the individual DOFs. Once the
   equation numbers have been set the numberer then invokes `setID()` on all the
   `FE_Elements` in the model. Finally the `DOF_Numberer` invokes `setNumEqn()`
   on `theAnalysisModel`.

4. Then `doneNumberingDOF()` is invoked on the `ConstraintHandler` which
   invokes `setID()` on all the `FE_Elements` in the model. 

5. It invokes `setSize(theModel.getDOFGraph())`{.cpp} on `theSOE` and
   `theEigenSOE` which causes the system of equation to determine its size
   based on the connectivity of the dofs in the analysis model.

6. It then invokes `domainChanged()` on `theIntegrator` and
   `theAlgorithm` to inform these objects that changes have occurred in
   the model.

Returns $0$ if successful. At any stage above, if an
error occurs the method is stopped, a warning message is printed and
a negative number is returned.


### Public Methods to vary the type of Analysis

:::{.admonition}
```cpp
int setAlgorithm(EquiSolnAlgo &theNewAlgorithm)
```
:::
To change the algorithm between analysis. It first invokes the
destructor on the old `SolutionAlgorithm` object associated with the
analysis. It then sets the `SolutionAlgorithm` associated with the
analysis to be `newAlgorithm` and sets the links for this object by
invoking `setLinks()`. Checks then to see if the domain has changed, if
true it invokes `domainChanged()`, otherwise it invokes
`domainChanged()` on the new SolutionAlgorithm. Returns $0$ if
successful, a warning message and a negative number if not.

:::{.admonition}
```cpp
int setIntegrator(TransientIntegrator &theNewIntegrator)
```
:::
To change the integration scheme between analysis. It first invokes the
destructor on the old Integrator object associated with the analysis. It
then sets the SolutionAlgorithm associated with the analysis to be
`newAlgorithm` and sets the links for this object by invoking
`setLinks()`. It also invokes `setLinks()` on the ConstraintHandler and
SolutionAlgorithm objects. Checks then to see if the domain has changed,
if true it invokes `domainChanged()`, otherwise it invokes
`domainChanged()` on the new Integrator. Returns $0$ if successful, a
warning message and a negative number if not.

:::{.admonition}
```cpp
int setLinearSOE(LinearSOE &theNewSOE);
int setEigenSOE(EigenSOE &theNewSOE);
```
:::
To change the respective system of equation object between analysis. It
first invokes the destructor on the old LinearSOE object associated with
the analysis. It then sets the SolutionAlgorithm associated with the
analysis to be `newSOE`. links for this object by invoking `setLinks()`.
It then invokes `setLinks()` on the ConstraintHandler and
SolutionAlgorithm objects. Checks then to see if the domain has changed,
if true it invokes `domainChanged()`, otherwise it invokes `setSize()`
on the new LinearSOE. Returns $0$ if successful, a warning message and a
negative number if not.

