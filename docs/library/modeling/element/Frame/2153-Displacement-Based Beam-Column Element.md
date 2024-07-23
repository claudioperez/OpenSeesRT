# Displacement-Based Beam-Column

This command is used to construct a displacement beam element object,
which is based on the displacement formulation, and considers the spread
of plasticity along the element.

```tcl
element dispBeamColumn $eleTag $iNode $jNode $numIntgrPts
        $secTag $transfTag < -mass $massDens > < -cMass >
        < -integration $intType >
```

<p>To change the sections along the element length, the following form
of command may be used:</p>

```tcl
element dispBeamColumn $eleTag $iNode $jNode $numIntgrPts
        -sections $secTag1 $secTag2 ... $transfTag < -mass $massDens >
        < -cMass > < -integration $intType >
```

<hr />
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
<td><code class="parameter-table-variable">numIntgrPts</code></td>
<td><p>number of integration points along the element.</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">secTag</code></td>
<td><p>identifier for previously-defined section object</p></td>
</tr>
<tr class="odd">
<td><p><strong>`secTag1` `secTag2` ...</strong></p></td>
<td><p>`numIntgrPts` identifiers of previously-defined section
object</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">transfTag</code></td>
<td><p>identifier for previously-defined coordinate-transformation
(CrdTransf) object</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">massDens</code></td>
<td><p>element mass density (per unit length), from which a lumped-mass
matrix is formed (optional, default = 0.0)</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-flag">-cMass</code></p></td>
<td><p>to form consistent mass matrix (optional, default = lumped mass
matrix)</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">intType</code></td>
<td><p>numerical integration type, options are Lobotto, Legendre, Radau,
NewtonCotes, Trapezoidal (optional, default = Legendre)</p></td>
</tr>
</tbody>
</table>

<p>NOTE:</p>
<ul>
<li>The default integration along the element is based on Gauss-Legendre quadrature rule.</li>

<li>The default element is prismatic, i.e. the beam is represented by
the section model identified by `secTag` at each integration point.</li>

<li>The valid queries to a displacement-based beam-column element when
  creating an ElementRecorder object are `force`, and 
  `section $secNum $secArg1 $secArg2...` Where `secNum` refers to the integration point whose
  data is to be output valid entries being 1 through `numIntgrPts`.</li>
</ul>

## Examples

displacement-based beam column
element added with tag 1 between nodes 2 and 4 that has 5 integration
points, each using section 8, and the element uses geometric
transformation `9`

```tcl
element dispBeamColumn 1 2 4 5 8 9; 
```

## References
<hr />
<p>Code Developed by: <span style="color:blue"> Michael H. Scott,
Oregon State University </span></p>

