# Elastic Timoshenko Beam Column

<p>This command is used to construct an ElasticTimoshenkoBeam element
object. A Timoshenko beam is a frame member that accounts for shear
deformations. The arguments for the construction of an elastic
Timoshenko beam element depend on the dimension of the problem, ndm:</p>
<p>For a two-dimensional problem:</p>

```tcl
element ElasticTimoshenkoBeam $eleTag $iNode $jNode $E $G
        $A $Iz $Avy $transfTag &lt;-mass $massDens&gt;
        &lt;-cMass&gt;
```

<p>For a three-dimensional problem:</p>

```tcl
element ElasticTimoshenkoBeam $eleTag $iNode $jNode $E $G
        $A $Jx $Iy $Iz $Avy $Avz $transfTag &lt;-mass $massDens&gt;
        &lt;-cMass&gt;
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
<td><code class="parameter-table-variable">E</code></td>
<td><p>Young's Modulus</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">G</code></td>
<td><p>Shear Modulus</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">A</code></td>
<td><p>cross-sectional area of element</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">Jx</code></td>
<td><p>torsional moment of inertia of cross section</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">Iy</code></td>
<td><p>second moment of area about the local y-axis</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">Iz</code></td>
<td><p>second moment of area about the local z-axis</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">Avy</code></td>
<td><p>Shear area for the local y-axis</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">Avz</code></td>
<td><p>Shear area for the local z-axis</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">transfTag</code></td>
<td><p>identifier for previously-defined coordinate-transformation
(CrdTransf) object</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">massDens</code></td>
<td><p>element mass per unit length (optional, default = 0.0)</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-flag">-cMass</code></p></td>
<td><p>to form consistent mass matrix (optional, default = lumped mass
matrix)</p></td>
</tr>
</tbody>
</table>
<p>NOTES:</p>
<p>The valid queries to an elastic Timoshenko beam element when creating
an ElementRecorder object are 'force'.</p>
<p>For solid rectangular sections, the shear area is 5/6 of the gross
area. For solid circular sections, the shear area is 9/10 of the gross
area. For I-shapes, the shear area can be approximated as Aweb.</p>

## Examples

<p>element ElasticTimoshenkoBeam 1 2 4 100.0 45.0 6.0 4.5 5.0 9; #
elastic Timoshenko element with tag 1 between nodes 2 and 4 with E=100,
G=45, A=6.0, I=4.5 and Av=5.0 which uses transformation 9</p>
<hr />
<p>Code Developed by: <span style="color:blue"> Andreas
Schellenberg </span></p>
