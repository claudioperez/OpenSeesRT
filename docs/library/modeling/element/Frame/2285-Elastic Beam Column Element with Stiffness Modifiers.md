---
descrption: Elastic Beam Column Element with Stiffness Modifiers
---
# ModElasticBeam2d

<p>This command is used to construct a ModElasticBeam2d element object.
The arguments for the construction of an elastic beam-column element
with stiffness modifiers is applicable for 2-D problems. This element
should be used for modelling of a structural element with an equivalent
combination of one elastic element with stiffness-proportional damping,
and two springs at its two ends with no stiffness proportional damping
to represent a prismatic section. The modelling technique is based on a
number of analytical studies discussed in Zareian and Medina (2010) and
Zareian and Krawinkler (2009) and is utilized in order to solve problems
related to numerical damping in dynamic analysis of frame structures
with concentrated plasticity springs.</p>

```tcl
element ModElasticBeam2d $eleTag $iNode $jNode $A $E $Iz
        $K11 $K33 $K44 $transfTag &lt;-mass $massDens&gt;
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
<td><code class="parameter-table-variable">A</code></td>
<td><p>cross-sectional area of element</p></td>
</tr>
<tr class="even">
<td><code class="parameter-table-variable">E</code></td>
<td><p>Young's Modulus</p></td>
</tr>
<tr class="odd">
<td><code class="parameter-table-variable">Iz</code></td>
<td><p>second moment of area about the local z-axis</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">K11</code></p></td>
<td><p>stiffness modifier for translation</p></td>
</tr>
<tr class="odd">
<td><p><code class="parameter-table-variable">K33</code></p></td>
<td><p>stiffness modifier for translation</p></td>
</tr>
<tr class="even">
<td><p><code class="parameter-table-variable">K44</code></p></td>
<td><p>stiffness modifier for rotation</p></td>
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
<p><strong>Element Formation</strong>:</p>
<hr />
<p>For structural elements with time invariant moment gradient, a
two-dimensional, prismatic beam element with six degrees of freedom is
to be replaced with a two-dimensional, prismatic beam element composed
of semi-rigid rotational springs at the ends and an elastic beam element
in the middle. The rotational stiffness at the end of the original beam
element is Ke = 6EIz/L (where E is the modulus of elasticity, Iz the
moment of inertia, and L the length of the beam), and the ratio of the
rotational spring stiffness, Ks, to the elastic beam stiffness, Ke, of
the modified beam element is defined as n = Ks/Ke.</p>
<hr />
<table>
<tbody>
<tr class="odd">
<td></td>
</tr>
</tbody>
</table>
<p><strong>Elastic Element with 2-end Springs:</strong></p>
<table>
<tbody>
<tr class="odd">
<td><p>1. The elastic element in between the two springs should have an
elastic moment of inertia equal to Iz,mod = (n+1)/n * Iz</p></td>
</tr>
<tr class="even">
<td><p>2. The "n" times stiff rotational springs should have an elastic
stiffness of Ks=n*6*EIz,mod/L</p></td>
</tr>
<tr class="odd">
<td><p>3. the elastic element should have an elastic stiffness
coefficient K44 = 6*(1+n)/(2+3*n)</p></td>
</tr>
<tr class="even">
<td><p>4. The elastic element should have an elastic stiffness
coefficient K11 = K33 = (1+2*n)*K44/(1+n)</p></td>
</tr>
<tr class="odd">
<td><p>5. The modified stiffness coefficient <em>bmod</em> for stiffness
proportional damping of the elastic element is: "bmod" =
1+(1/2n)*b</p></td>
</tr>
<tr class="even">
<td></td>
</tr>
</tbody>
</table>
<hr />
<table>
<tbody>
<tr class="odd">
<td></td>
</tr>
</tbody>
</table>
<p><strong>Elastic Element with 1-end Spring:</strong></p>
<table>
<tbody>
<tr class="odd">
<td><p>1. The elastic element in between the two springs should have an
elastic moment of inertia equal to Iz,mod = (n+1)/n * Iz</p></td>
</tr>
<tr class="even">
<td><p>2. The "n" times stiff rotational springs should have an elastic
stiffness of Ks=n*6*EIz,mod/L</p></td>
</tr>
<tr class="odd">
<td><p>3. The elastic element should have an elastic stiffness
coefficient K44 = 6*n/( 1+3*n )</p></td>
</tr>
<tr class="even">
<td><p>4. the elastic element should have an elastic stiffness
coefficient K11 = (1+2*n)*K44/(1+n)</p></td>
</tr>
<tr class="odd">
<td><p>5. the elastic element should have an elastic stiffness
coefficient K33 = 2*K44</p></td>
</tr>
<tr class="even">
<td><p>6. The modified stiffness coefficient <em>bmod</em> for stiffness
proportional damping of the elastic element is: "bmod" =
1+(1/2n)*b</p></td>
</tr>
<tr class="odd">
<td></td>
</tr>
</tbody>
</table>
<hr />
<table>
<tbody>
<tr class="odd">
<td></td>
</tr>
</tbody>
</table>
<p><strong>EXAMPLE</strong>:</p>
<table>
<tbody>
<tr class="odd">
<td><p>element ModelasticBeam2d 1 2 4 5.5 100.0 1e6 4.0 4.0 2.0 1; #
elastic element tag 1 between nodes 2 and 4 with area 5.5, E 100 and IZ
1e6 with K11=K33=4.0, K44=2.0 which uses transformation 1</p></td>
</tr>
<tr class="even">
<td><p><strong>Note:</strong> For n = 1 rigid spring, the stiffness
coefficients of the elastic springs are as follows: K11=K33 = 4.0 and
K44 = 2.0</p></td>
</tr>
<tr class="odd">
<td></td>
</tr>
</tbody>
</table>
<hr />
<p><strong>References</strong>:</p>
<table>
<tbody>
<tr class="odd">
<td><p><strong>[1]</strong></p></td>
<td><p>Zareian, F. and Medina, R. A. (2010). “A practical method for
proper modeling of structural damping in inelastic plane structural
systems,” Computers &amp; Structures, Vol. 88, 1-2, pp. 45-53.</p></td>
</tr>
<tr class="even">
<td><p><strong>[2]</strong></p></td>
<td><p>Zareian, F. and Krawinkler, H. (2009). "Simplified
performance-based earthquake engineering" Technical Report 169, The John
A. Blume Earthquake Engineering Research Center, Department of Civil
Engineering, Stanford University, Stanford, CA. [electronic version: <a
href="https://blume.stanford.edu/tech_reports">https://blume.stanford.edu/tech_reports</a>]</p></td>
</tr>
<tr class="odd">
<td></td>
<td></td>
</tr>
</tbody>
</table>
<hr />
<p>Code Developed by: <span style="color:blue"> by Dr. Dimitrios
G. Lignos (McGill University) </span></p>
