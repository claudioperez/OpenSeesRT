# Force Formulation

The ForceBeamColumn element is based on the iterative force-based formulation. 
A variety of numerical integration options can be used in the element state
determination and encompass both distributed plasticity and plastic
hinge integration. See <a href="image:IntegrationTypes.pdf"
title="wikilink">image:IntegrationTypes.pdf</a> for more details on the
available numerical integration options.

:::{apidoc="opensees.element.ForceBeamColumn"}
```tcl
element ForceBeamColumn $eleTag $iNode $jNode $transfTag
  { IntegrationType arg1 arg2 ... }
  < -mass $massDens > < -iter  $maxIters $tol >
```

<table>
<tbody>
<tr class="odd">
<td><code class="parameter-table-variable">eleTag</code></td>
<td><p>unique element object tag</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">iNode jNode</code></p></td>
<td><p>end nodes</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">transfTag</code></td>
<td><p>identifier for previously-defined coordinate-transformation
(CrdTransf) object</p></td>
</tr>
<tr class="even">
<td><p><strong>IntegrationType arg1 arg2 ...</strong></p></td>
<td><p>specifies locations and weights of integration points and their
associated section force-deformation models (see <a
href="image:IntegrationTypes.pdf"
title="wikilink">image:IntegrationTypes.pdf</a>)</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">massDens</code></td>
<td><p>element mass density (per unit length), from which a lumped-mass
matrix is formed (optional, default=0.0)</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">maxIters</code></td>
<td><p>maximum number of iterations to undertake to satisfy element
compatibility (optional, default=10)</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">tol</code></td>
<td><p>tolerance for satisfaction of element compatibility (optional,
default=10-12)</p></td>
</tr>
</tbody>
</table>
:::

<p>Original command that assumes Gauss-Lobatto integration with a copy
of the same section force-deformation model at each integration
point:</p>

```tcl
element forceBeamColumn $eleTag $iNode $jNode
        $numIntgrPts $secTag $transfTag < -mass $massDens > < -iter
        $maxIters $tol > < -integration $intType >
```

<table>
<tbody>
<tr class="odd">
<td><code class="parameter-table-variable">eleTag</code></td>
<td><p>unique element object tag</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">numIntgrPts</code></td>
<td><p>number of Gauss-Lobatto integration points along the
element.</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">secTag</code></td>
<td><p>identifier for previously-defined section object</p></td>
</tr>
</tbody>
</table>
<p>Alternative command (kept for backward compatability):</p>

```tcl
element nonlinearBeamColumn $eleTag $iNode $jNode
        $numIntgrPts $secTag $transfTag < -mass $massDens > 
        < -iter $maxIters $tol > < -integration $intType >
```

<table>
<tbody>
<tr class="odd">
<td><code class="parameter-table-variable">eleTag</code></td>
<td><p>unique element object tag</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">intType</code></td>
<td><p>numerical integration type, options are Lobatto, Legendre, Radau,
NewtonCotes, Trapezoidal (optional, default= Lobatto)</p></td>
</tr>
</tbody>
</table>
<hr />

<p>NOTE:</p>
<p>The following three commands give the same element definition (with
Gauss-Lobatto integration) despite some apparent permutations of the
input arguments:</p>
<ol>
<li>element forceBeamColumn $eleTag $iNode $jNode $transfTag Lobatto
$secTag $numIntgrPts</li>
<li>element forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag
$transfTag</li>
<li>element nonlinearBeamColumn $eleTag $iNode $jNode $numIntgrPts
$secTag $transfTag</li>
</ol>

## Notes
- The `-iter` switch enables the iterative form of the flexibility
  formulation. Note that the iterative form can improve the rate of global
  convergence at the expense of more local element computation.
- The valid response elements that an element of this type will
  respond to are:

  - `force` or `globalForce`
  - `localForce`
  - `basicForce`
  - `section $sectionNumber $arg1 $arg2 ...` (note: `sectionNumer` is
    integer 1 through `numIntegrPts`)
  - `basicDeformation`
  - `plasticDeformation`
  - `inflectionPoint`
  - `tangentDrift`
  - `integrationPoints`
  - `integrationWeights`

- Here is a link to the source code to obtain information about the
  location and weight of the Gauss-Lobatto integration points <a
  href="http://opensees.berkeley.edu/WebSVN/filedetails.php?repname=OpenSees&path=%2Ftrunk%2FSRC%2Felement%2FforceBeamColumn%2FLobattoBeamIntegration.cpp">1</a>


## Examples

force beam column element added with tag 1 between nodes 2 and 4 that has Gauss-Lobatto 5
integration points, each using section 8, and the element uses geometric
transformation 9
```tcl
element forceBeamColumn 1 2 4 9 Lobatto 8 5; 
```

<p>FURTHER DOCUMENTATION ON INTEGRATION OPTIONS:</p>
<p><a href="image:IntegrationTypes.pdf"
title="wikilink">image:IntegrationTypes.pdf</a></p>
## References
<ul>
<li>Neuenhofer, Ansgar, FC Filippou. Geometrically Nonlinear
Flexibility-Based Frame Finite Element. ASCE Journal of Structural
Engineering, Vol. 124, No. 6, June, 1998. ISSN
0733-9445/98/0006-0704-0711. Paper 16537. pp. 704-711.</li>
</ul>
<ul>
<li>Neuenhofer, Ansgar, FC Filippou. Evaluation of Nonlinear Frame
Finite-Element Models. ASCE Journal of Structural Engineering, Vol. 123,
No. 7, July, 1997. ISSN 0733-9445/97/0007-0958-0966. Paper No. 14157.
pp. 958-966.</li>
</ul>
<ul>
<li>Neuenhofer, Ansgar, FC Filippou. ERRATA -- Geometrically Nonlinear
Flexibility-Based Frame Finite Element. ASCE Journal of Structural
Engineering, Vol. 124, No. 6, June, 1998. ISSN
0733-9445/98/0006-0704-0711. Paper 16537. pp. 704-711.</li>
</ul>
<ul>
<li>Taucer, Fabio F, E Spacone, FC Filippou. A Fiber Beam-Column Element
for Seismic Response Analysis of Reinforced Concrete Structures. Report
No. UCB/EERC-91/17. Earthquake Engineering Research Center, College of
Engineering, University of California, Berkeley. December 1991.</li>
</ul>
<ul>
<li>Spacone, Enrico, V Ciampi, FC Filippou. A Beam Element for Seismic
Damage Analysis. Report No. UCB/EERC-92/07. Earthquake Engineering
Research Center, College of Engineering, University of California,
Berkeley. August 1992.</li>
</ul>
<hr />

Code maintained by: <a>cmp</a>

