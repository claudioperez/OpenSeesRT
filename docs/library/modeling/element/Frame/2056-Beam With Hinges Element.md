# BeamWithHinges

This command is used to construct a forceBeamColumn element object,
which is based on the non-iterative (or iterative) flexibility
formulation. The locations and weights of the element integration points
are based on so-called plastic hinge integration, which allows the user
to specify plastic hinge lenghts at the element ends. Two-point Gauss
integration is used on the element interior while two-point Gauss-Radau
integration is applied over lengths of $4\ell_{p,i}$ and $4\ell_{p,j}$ at the element
ends, viz. "modified Gauss-Radau plastic hinge integration". A total of
six integration points are used in the element state determination (two
for each hinge and two for the interior).

<p>Users may be familiar with the <strong>beamWithHinges</strong>
command format (see below); however, the format shown here allows for
the simple but important case of using a material nonlinear section
model on the element interior. The previous beamWithHinges command
constrained the user to an elastic interior, which often led to
unconservative estimates of the element resisting force when plasticity
spread beyond the plastic hinge regions in to the element interior.</p>

The advantages of this new format over the previous
<strong>beamWithHinges</strong> command are

<ul>
  <li>Plasticity can spread beyond the plastic hinge regions</li>
  <li>Hinges can form on the element interior, e.g., due to distributed member loads</li>
</ul>

```tcl
element forceBeamColumn $eleTag $iNode $jNode $transfTag
        "HingeRadau $secTagI $LpI $secTagJ $LpJ $secTagInterior" 
        < -mass $massDens > 
        < -iter $maxIters $tol >
```

<table>
<tbody>
<tr class="odd">
<td><code class="parameter-table-variable">eleTag</code></td>
<td><p>unique element object tag</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">iNode jNode</code></p></td>
<td><p>nodes at element ends I and J, respectively</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">secTagI</code></td>
<td><p>identifier for previously-defined section object at end I</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">LpI</code></td>
<td><p>plastic hinge length at end I</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">secTagJ</code></td>
<td><p>identifier for previously-defined section object at end J</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">LpJ</code></td>
<td><p>plastic hinge length at end J</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">secTagInterior</code></td>
<td><p>identifier for previously-defined section object on the element
interior (DOES NOT HAVE TO BE ELASTIC, but can be any type of section,
including fiber)</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">transfTag</code></td>
<td><p>identifier for previously-defined
coordinate-transformation</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">maxIters</code></td>
<td><p>maximum number of iterations to undertake to satisfy element
compatibility (optional, default=1)</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">tol</code></td>
<td><p>tolerance for satisfaction of element compatibility (optional,
default=10-16)</p></td>
</tr>
</tbody>
</table>

NOTE: The keyword <strong>HingeRadau</strong> can be changed to one
of the following in order to use a different hinge integration
approach:

<ul>
  <li><strong>HingeRadau</strong> -- two-point Gauss-Radau applied to the hinge regions over $4\ell_{pI}$ and $4\ell_{pJ}$ (six element integration points)</li>
  <li><strong>HingeRadauTwo</strong> -- two-point Gauss-Radau in the hinge regions applied over $\ell_{pI}$ and $\ell_{pJ}$ (six element integration points)</li>
  <li><strong>HingeMidpoint</strong> -- midpoint integration over the hinge regions (four element integration points)</li>
  <li><strong>HingeEndpoint</strong> -- endpoint integration over the hinge regions (four element integration points)</li>
</ul>

For more information on the behavior, advantages, and disadvantages
of these approaches to plastic hinge integration, see

- Scott, M.H. and G.L. Fenves. "<a
  href="http://dx.doi.org/10.1061/(ASCE)0733-9445(2006)132:2(244)">Plastic
  Hinge Integration Methods for Force-Based Beam-Column Elements</a>",
  Journal of Structural Engineering, 132(2):244-252, February 2006.

- Scott, M.H. and K.L. Ryan. "<a
  href="http://dx.doi.org/10.1193/1.4000136">Moment-Rotation Behavior of
  Force-Based Plastic Hinge Elements</a>", Earthquake Spectra, 29(2):597-607, May 2013.

<p>The primary advantages of <strong>HingeRadau</strong> are</p>
<ul>
  <li>The user can specify a physically meaningful plastic hinge length</li>
  <li>The largest bending moment is captured at the element ends</li>
  <li>The exact numerical solution is recovered for a linear-elastic prismatic beam</li>
  <li>The characteristic length is equal to the user-specified plastic hinge length when deformations localize at the element ends</li>
</ul>

<p>while the primary disadvantages are</p>
<ul>
  <li>The element post-yield response is too flexible for strain-hardening section response (consider using <strong>HingeRadauTwo</strong>)</li>
  <li>The user needs to know the plastic hinge length <em>a priori</em> (empirical equations are available)</li>
</ul>

<p>NOTE: See the <strong><a href="Force-Based_Beam-Column_Element"
title="wikilink"> forceBeamColumn</a></strong> page for valid recorder queries.</p>
<hr />

<p>Original command (maintained for backward compatibility)</p>
<p>NOTE: this form of the command forces the element interior to be
linear-elastic, which is not always the best approach.</p>
<p>For 2D:</p>

```tcl
element beamWithHinges $eleTag $iNode $jNode $secTagI
        $Lpi $secTagJ $Lpj $E $A $Iz $transfTag < -mass $massDens >
        < -iter $maxIters $tol >
```

<p>For 3D:</p>

```tcl
element beamWithHinges $eleTag $iNode $jNode $secTagI
        $Lpi $secTagJ $Lpj $E $A $Iz $Iy $G $J $transfTag < -mass
        $massDens > < -iter $maxIters $tol >
```

<p>All inputs are the same as above, with the following additional
inputs, which are used solely to create a "dummy" elastic section at the
two Gauss integration points of the element interior</p>
<table>
<tbody>
<tr class="odd">
<td><code class="parameter-table-variable">E</code></td>
<td><p>Youngs modulus elastic portion</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">A</code></td>
<td><p>Area for elastic portion</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">Iz</code></td>
<td><p>second moment of area for elastic portion about local z</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">Iy</code></td>
<td><p>second moment of area for elastic portion about local y</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">G</code></td>
<td><p>torsional moment of inertia of cross section for elastic
portion</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">J</code></td>
<td><p>Shear Modulus of elastic portion.</p></td>
</tr>
</tbody>
</table>

## C++ Interface (BeamWithHinges2D)

```cpp
#include <element/BeamWithHinges2D.h>

class BeamWithHinges2D: public DomainComponent
```


BeamWithHinges2D is a beam-column element which uses the force based
formulation for its state determination. This element has material
non-linear hinges at both ends and exhibits linear elastic behavior
through its interior region, including linear elastic shear effects.

### Constructors

### Destructor


Constructs a BeamWithHinges2D object with tag `tag`, unique among all
elements in the domain. The end nodes of the element are set to be those
with tags `nodeI` and `nodeJ`. The elastic properties of the interior
beam region are set to be `E`, the modulus of elasticity; `I`, the
second moment of area of the beam cross-section; `A`, the beam
cross-sectional area; `G`, the shear modulus; and `alpha`, the
cross-section shape factor for elastic shear effects. Element hinges are
created by obtaining copies of `sectionI` and `sectionJ`. The hinge
lengths are specified as ratios of the total element length, `ratioI`
and `ratioJ`. The element distributed load (reference value) is set to
be `distrLoad`, and the element mass density per unit length is
`massDens`.

Constructs a BeamWithHinges2D object with tag `tag`, unique among all
elements in the domain. The end nodes of the element are set to be those
with tags `nodeI` and `nodeJ`. The elastic properties of the interior
beam region are set to be `E`, the modulus of elasticity; `I`, the
second moment of area of the beam cross-section; `A`, the beam
cross-sectional area; `G`, the shear modulus; and `alpha`, the
cross-section shape factor for elastic shear effects. Element hinges are
created by obtaining copies of `sectionI` and `sectionJ`. The hinge
lengths are specified as ratios of the total element length, `ratioI`
and `ratioJ`. The element mass density per unit length is `massDens`.



Constructs a blank BeamWithHinges2D object.


Invokes the section destructors.


// Public Methods dealing with Nodes and dof

\

\
// Public Methods dealing with State

\

\
// Public Methods for obtaining Linearized Stiffness, Mass and Damping Matrices

\

\

// Public Methods for obtaining Resisting Forces\

\

\

\
Returns 2, the number of external nodes for this element.

Returns an ID containing the tags of the two external nodes for this
element.

Returns 6, the number of degrees of freedom for this element.

```{.cpp}
int commitState(void);
```

Invokes `commitState()` on each section and returns the sum of the
result of these invocations.

```{.cpp}
int revertToLastCommit(void);
```

Invokes `revertToLastCommit()` on each section and returns the sum of
the result of these invocations.

```{.cpp}
int revertToStart(void);
```

Invokes `revertToStart()` on each section and returns the sum of the
result of these invocations.

```{.cpp}
const Matrix &getTangentStiff(void);
```

Computes the element flexibility matrix, then returns its inverse, the
element stiffness matrix. The element flexibility is the sum of the
hinge flexibilities, $\mathbf{f}_I$ and $\mathbf{f}_J$, and the elastic
flexibility of the beam interior, $\mathbf{f}_{mid}$.

$$\label{eq:fele}
\fbas := \int_{0}^{L}{\bint^T \fsec \bint \: dx} = \mathbf{f}_I + \mathbf{f}_{mid} + \mathbf{f}_J$$

The flexibility of the beam interior is obtained in closed form,

$$\mathbf{f}_{mid} = \int_{l_I}^{L-l_J}{\bint^T \fsec_{mid} \bint \: dx}$$

where $\bint$ is the force interpolation matrix,

$$\bint := \left[
   \begin{array}{c c c}
      1 &           0 &               0 \\
      0 & \frac{x}{L} & \frac{x}{L} - 1 \\
      0 & \frac{1}{L} &     \frac{1}{L}
   \end{array} 
 \right]$$

and $\fsec$ is the elastic section flexibility of the beam interior.

$$\fsec_{mid} = \left[
   \begin{array}{c c c}
      \frac{1}{EA} &            0 &                   0 \\
                 0 & \frac{1}{EI} &                   0 \\
                 0 &            0 & \frac{1}{\alpha GA}
   \end{array}
 \right]$$

The hinge flexibilities, $\mathbf{f}_I$ and $\mathbf{f}_J$, are obtained
by the midpoint integration rule,

$${\mathbf{f}}_i = \bint(x_i)^T \fsec_i \bint(x_i) * l_i, \:\: i=I,J$$

where $x_i$ is the midpoint of hinge $i$, measured along the length of
the beam, and is the point at which the force interpolation matrix,
$\bint$ is evaluated. The flexiblity, $\fsec_i$, is obtained from the
constitutive relation for section $i$.
The element stiffness is then obtained by inversion of the element
flexibility, given by Equation
[\[eq:fele$$
](#eq:fele){reference-type="ref" reference="eq:fele"}.

$$\label{eq:kele}
\kbas = \fbas^{-1}$$

The element then obtains the matrix, $\mathbf{A}$, which transforms the
element basic stiffness from its corotating frame to the global frame of
reference. The transformed stiffness matrix, $\kele$, is then assembled
into the structural system of equations.

$$\kele = \mathbf{A}^T \kbas \mathbf{A}$$\

```{.cpp}
const Matrix &getSecantStiff(void);
```

To return the elements secant stiffness matrix. THIS SECANT MAY BE
REMOVED.

```{.cpp}
const Matrix &getDamp(void);
```

To return the damping matrix.

```{.cpp}
const Matrix &getMass(void);
```

Returns the element lumped mass matrix, $\textrm{mass}$. It is assumed that the
mass density per unit length, $\rho$, is constant along the entire
element, including the hinge regions.

$$\label{eq:mele}
\textrm{mass} = \left[
   \begin{array}{c c c c c c}
      \frac{\rho L}{2} & 0 & 0 & 0 & 0 & 0 \\
      0 & \frac{\rho L}{2} & 0 & 0 & 0 & 0 \\
      0 & 0 & 0 & 0 & 0 & 0 \\
      0 & 0 & 0 & \frac{\rho L}{2} & 0 & 0 \\
      0 & 0 & 0 & 0 & \frac{\rho L}{2} & 0 \\
      0 & 0 & 0 & 0 & 0 & 0
   \end{array}
 \right]$$


```{.cpp}
void zeroLoad(void);
```

This is a method invoked to zero the element load contributions to the
residual.

```{.cpp}
const Vector &getResistingForce(void);
```

Returns the element resisting force vector. The basic element force
vector is obtained as the product of the basic element stiffness,
$\kbas$, given by Equation [\[eq:kele$$
](#eq:kele){reference-type="ref"
reference="eq:kele"}, and the basic element deformations, $\vbas$.

$$\qbas = \kbas \vbas$$

The basic element force vector is then transformed from the corotating
frame to the global frame of reference. The transformed force vector is
then assembled into the structural system of equations.

$$\label{eq:qele}
\qele = \mathbf{A}^T \qbas$$\

```{.cpp}
const Vector &getResistingForceIncInertia(void);
```

Returns the element resisting force vector, $\tilde{\qele}$ with inertia
forces included,

$$\tilde{\qele} = \qele - \textrm{mass} \ddot{\mathbf u}$$

where $\qele$ and $\textrm{mass}$ are obtained from Equations
[\[eq:qele$$
](#eq:qele){reference-type="ref" reference="eq:qele"} and
[\[eq:mele$$
](#eq:mele){reference-type="ref" reference="eq:mele"},
respectively, and $\ddot{\mathbf u}$ is the vector of trial nodal
accelerations for the element.

<hr />

<p>Code maintained by: <a
href="http://web.engr.oregonstate.edu/~mhscott">Michael H. Scott, Oregon State University</a></p>

NEED TO ADD ADD_INERTIA_LOAD TO INTERFACE .. SEE EARTHQUAKE_PATTERN
CLASS.

