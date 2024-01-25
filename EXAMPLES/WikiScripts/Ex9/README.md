---
title:  Rectilinear Cross Sections
...

# Rectilinear Cross Sections

<div class="mw-parser-output">

<!-- <div id="toc" class="toc"> -->
<!-- <input type="checkbox" role="button" id="toctogglecheckbox" class="toctogglecheckbox" style="display:none"> -->
<!-- <div class="toctitle" lang="en" dir="ltr"> -->

For the case of the uniaxial section, moment-curvature and axial force-deformation curves are defined independently, and numerically.

For the case of the fiber sections (steel and RC), uniaxial materials are defined numerically (stress-strain relationship) and are combined into a fiber section where moment-curvature and axial force-deformation characteristics and their interaction are calculated computationally.

#### 2D vs. 3D

While this distinction does not affect the section definition itself, it
affects the degree-of-freedom associated with moment and curvature in the
subsequent analysis. There are two differences between the two models:

1. The space defined with the model command (Defining the model builder, `ndm=#dimension ndf=#dofs`)
2. In the 3D model, torsional stiffness needs to be aggregated to the section.


<table style="margin:0; background:none; border:3px solid #ccc;">
<tbody><tr>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<h4><span class="mw-headline" id="Uniaxial_Section">Uniaxial Section</span></h4>
<table style="width:100%; vertical-align:top;background:#white;">
<tbody><tr><td style="color:#000;">
<p><a href="/wiki/index.php/OpenSees_Example_9._Build_%26_Analyze_a_Section_Example" title="OpenSees Example 9. Build &amp; Analyze a Section Example"><img alt="ExampleSection Uniaxial.gif" src="./ExampleSection_Uniaxial.gif" width="306" height="215"></a>
</p>
</td></tr></tbody></table>
</td>
<!--
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Files</strong><ul><li><a href="./Ex9a.build.UniaxialSection2D.tcl" class="internal" title="Ex9a.build.UniaxialSection2D.tcl">Ex9a.build.UniaxialSection2D.tcl</a></li>
<li><a href="./Ex9a.build.UniaxialSection3D.tcl" class="internal" title="Ex9a.build.UniaxialSection3D.tcl">Ex9a.build.UniaxialSection3D.tcl</a></li>
<li><a href="./LibUnits.tcl" class="internal" title="LibUnits.tcl">LibUnits.tcl</a></li></ul>
</td>
-->

<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Notes</strong><ul><li>Flexure and axial behavior are uncoupled in this type of section</li></ul>
</td></tr></tbody></table>

<table style="margin:0; background:none; border:3px solid #ccc;">
<tbody><tr>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<h4><span class="mw-headline" id="Fiber_Section:_AISC_Standard_W_Section">Fiber Section: AISC Standard W Section</span></h4>
<table style="width:100%; vertical-align:top;background:#white;">
 <tbody><tr>
 <td style="color:#000;">
 <p><a href="/wiki/index.php/OpenSees_Example_9._Build_%26_Analyze_a_Section_Example" title="OpenSees Example 9. Build &amp; Analyze a Section Example"><img alt="ExampleSection W.gif" src="./ExampleSection_W.gif" width="312" height="163"></a>
 </p>
 </td></tr></tbody>
</table>
</td>

<!--
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Files</strong><ul><li><a href="./Ex9b.build.WSection2D.tcl" class="internal" title="Ex9b.build.WSection2D.tcl">Ex9b.build.WSection2D.tcl</a></li>
<li><a href="./Ex9b.build.WSection3D.tcl" class="internal" title="Ex9b.build.WSection3D.tcl">Ex9b.build.WSection3D.tcl</a></li>
<li><a href="./LibUnits.tcl" class="internal" title="LibUnits.tcl">LibUnits.tcl</a></li></ul>
</td>
-->
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Notes</strong><ul><li>Coupled biaxial flexure and axial behavior</li></ul>
</td></tr></tbody>
</table>

<table style="margin:0; background:none; border:3px solid #ccc;">
<tbody><tr>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<h4><span id="Fiber_Section:_Reinforced_Concrete_Section_--_Rectangular_Symmetric_Section,_Unconfined_Concrete"></span><span class="mw-headline" id="Fiber_Section:_Reinforced_Concrete_Section_--_Rectangular_Symmetric_Section.2C_Unconfined_Concrete">Fiber Section: Reinforced Concrete Section -- Rectangular Symmetric Section, Unconfined Concrete</span></h4><table style="width:100%; vertical-align:top;background:#white;">

<tbody><tr>
<td style="color:#000;">
<p><a href="/wiki/index.php/OpenSees_Example_9._Build_%26_Analyze_a_Section_Example" title="OpenSees Example 9. Build &amp; Analyze a Section Example"><img alt="ExampleSection RCSymmUnconf.gif" src="./ExampleSection_RCSymmUnconf.gif" width="419" height="221"></a>
</p>
</td></tr></tbody></table>
</td>
<!--
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Files</strong><ul><li><a href="./Ex9c.build.RCSection.RectUnconfinedSymm2D.tcl" class="internal" title="Ex9c.build.RCSection.RectUnconfinedSymm2D.tcl">Ex9c.build.RCSection.RectUnconfinedSymm2D.tcl</a></li>
<li><a href="./Ex9c.build.RCSection.RectUnconfinedSymm3D.tcl" class="internal" title="Ex9c.build.RCSection.RectUnconfinedSymm3D.tcl">Ex9c.build.RCSection.RectUnconfinedSymm3D.tcl</a></li>
<li><a href="./LibUnits.tcl" class="internal" title="LibUnits.tcl">LibUnits.tcl</a></li></ul>
</td>
-->

<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Notes</strong><ul><li>Coupled biaxial flexure and axial behavior</li></ul>
</td></tr></tbody>
</table>
<table style="margin:0; background:none; border:3px solid #ccc;">
<tbody><tr>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<h4><span id="Fiber_Section:_Reinforced_Concrete_Section_--_Rectangular_Symmetric_Section,_Confined_Concrete_Core"></span><span class="mw-headline" id="Fiber_Section:_Reinforced_Concrete_Section_--_Rectangular_Symmetric_Section.2C_Confined_Concrete_Core">Fiber Section: Reinforced Concrete Section -- Rectangular Symmetric Section, Confined Concrete Core</span></h4>
<table style="width:100%; vertical-align:top;background:#white;">
 <tbody><tr>
 <td style="color:#000;">
 <p><a href="/wiki/index.php/OpenSees_Example_9._Build_%26_Analyze_a_Section_Example" title="OpenSees Example 9. Build &amp; Analyze a Section Example"><img alt="ExampleSection RCSymmConf.gif" src="./ExampleSection_RCSymmConf.gif" width="476" height="223"></a>
 </p>
 </td></tr></tbody>
</table>
</td>
<!--
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Files</strong><ul><li><a href="./Ex9d.build.RCSection.RectConfinedSymm2D.tcl" class="internal" title="Ex9d.build.RCSection.RectConfinedSymm2D.tcl">Ex9d.build.RCSection.RectConfinedSymm2D.tcl</a></li>
<li><a href="./Ex9d.build.RCSection.RectConfinedSymm3D.tcl" class="internal" title="Ex9d.build.RCSection.RectConfinedSymm3D.tcl">Ex9d.build.RCSection.RectConfinedSymm3D.tcl</a></li>
<li><a href="./LibUnits.tcl" class="internal" title="LibUnits.tcl">LibUnits.tcl</a></li></ul>
</td>
-->
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Notes</strong><ul><li>Coupled biaxial flexure and axial behavior</li></ul>
</td></tr></tbody></table>
<table style="margin:0; background:none; border:3px solid #ccc;">
<tbody><tr>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<h4><span class="mw-headline" id="Fiber_Section:_Reinforced_Concrete_Section_--_Rectangular_Section">Fiber Section: Reinforced Concrete Section -- Rectangular Section</span></h4><table style="width:100%; vertical-align:top;background:#white;">

<tbody><tr>
<td style="color:#000;">
<p><a href="/wiki/index.php/OpenSees_Example_9._Build_%26_Analyze_a_Section_Example" title="OpenSees Example 9. Build &amp; Analyze a Section Example"><img alt="ExampleSection RCrect.gif" src="./ExampleSection_RCrect.gif" width="476" height="217"></a>
</p>
</td></tr></tbody></table>
</td>
<!--
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Files</strong><ul><li><a href="./Ex9e.build.RCSection.Rect2D.tcl" class="internal" title="Ex9e.build.RCSection.Rect2D.tcl">Ex9e.build.RCSection.Rect2D.tcl</a></li>
<li><a href="./Ex9e.build.RCSection.Rect3D.tcl" class="internal" title="Ex9e.build.RCSection.Rect3D.tcl">Ex9e.build.RCSection.Rect3D.tcl</a></li>
<li><a href="./LibUnits.tcl" class="internal" title="LibUnits.tcl">LibUnits.tcl</a></li></ul>
</td>
-->

<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Notes</strong><ul><li>Coupled biaxial flexure and axial behavior</li>
<li>generic rectangular section</li></ul>
<!--
<table style="width:100%; vertical-align:top;background:#white;">
<tbody><tr><td></td></tr></tbody></table>
-->
</td></tr></tbody></table>
<table style="margin:0; background:none; border:3px solid #ccc;">
<tbody><tr>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<h4><span id="Fiber_Section:_Reinforced_Concrete_Section_--_Circular_Section,_Confined_Core"></span><span class="mw-headline" id="Fiber_Section:_Reinforced_Concrete_Section_--_Circular_Section.2C_Confined_Core">Fiber Section: Reinforced Concrete Section -- Circular Section, Confined Core</span></h4>
<table style="width:100%; vertical-align:top;background:#white;">
<tbody><tr>
<td style="color:#000;">
<p><a href="/wiki/index.php/OpenSees_Example_9._Build_%26_Analyze_a_Section_Example" title="OpenSees Example 9. Build &amp; Analyze a Section Example"><img alt="ExampleSection RCcirc.gif" src="./ExampleSection_RCcirc.gif" width="459" height="167"></a>
</p>
</td></tr></tbody></table>
</td>
<!--
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Files</strong><ul><li><a href="./Ex9f.build.RCSection.Circ2D.tcl" class="internal" title="Ex9f.build.RCSection.Circ2D.tcl">Ex9f.build.RCSection.Circ2D.tcl</a></li>
<li><a href="./Ex9f.build.RCSection.Circ3D.tcl" class="internal" title="Ex9f.build.RCSection.Circ3D.tcl">Ex9f.build.RCSection.Circ3D.tcl</a></li>
<li><a href="./LibUnits.tcl" class="internal" title="LibUnits.tcl">LibUnits.tcl</a></li></ul>
</td>
-->
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Notes</strong><ul><li>Coupled biaxial flexure and axial behavior</li>
<li>generic circular section</li></ul>
<!--
<table style="width:100%; vertical-align:top;background:#white;">
<tbody><tr><td></td></tr></tbody></table>
-->
</td></tr></tbody></table>
<table style="margin:0; background:none; border:3px solid #ccc;">
<tbody><tr>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<h4><span id="Fiber_Section:_Reinforced_Concrete_Hollow_Section_--_Symmetric_Section,_Confined_Concrete"></span><span class="mw-headline" id="Fiber_Section:_Reinforced_Concrete_Hollow_Section_--_Symmetric_Section.2C_Confined_Concrete">Fiber Section: Reinforced Concrete Hollow Section -- Symmetric Section, Confined Concrete</span></h4>
<table style="width:100%; vertical-align:top;background:#white;">
  <tbody><tr><td style="color:#000;">
  <p><a href="/wiki/index.php/OpenSees_Example_9._Build_%26_Analyze_a_Section_Example" title="OpenSees Example 9. Build &amp; Analyze a Section Example"><img alt="HollowRC.png" src="./HollowRC.png" width="302" height="349"></a></p>
  </td></tr></tbody>
</table>
</td>

<!--
<td style="margin:0; width:25%; background:#white; vertical-align:top;">

<strong>Files</strong><ul><li><a href="./Ex9g.build.HollowSection3D.tcl" class="internal" title="Ex9g.build.HollowSection3D.tcl">Ex9g.build.HollowSection3D.tcl</a></li>
<li><a href="./LibUnits.tcl" class="internal" title="LibUnits.tcl">LibUnits.tcl</a></li></ul>

<table style="width:100%; vertical-align:top;background:#white;">


<tbody><tr><td></td></tr></tbody></table>
</td>

-->

<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Notes</strong><ul><li>Coupled biaxial flexure and axial behavior</li></ul>
<!--
<table style="width:100%; vertical-align:top;background:#white;">
<tbody><tr><td></td></tr></tbody></table>
-->
</td></tr></tbody></table>
<h3><span class="mw-headline" id="Moment-Curvature_Analysis">Moment-Curvature Analysis</span></h3>
<ul><li>This example introduces the moment-curvature procedures for sections in 2D or 3D space, as built in the previous section. (the only difference between them is the degree-of-freedom corresponding to curvature).</li>
<li>The moment-curvature analysis of a section is by creating a zero-length rotational-spring element. This section is subjected to a user-defined constant axial load and to a linearly-increasing moment to a user-defined maximum curvature.</li></ul>
<!-- <p><br></p> -->

<table style="margin:0; background:none; border:3px solid #ccc;">
<tbody><tr>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<h4><span class="mw-headline" id="2D_Moment-Curvature_Analysis">2D Moment-Curvature Analysis</span></h4>
<!--
<table style="width:100%; vertical-align:top;background:#white;">
<tbody><tr><td style="color:#000;"></td></tr></tbody>
</table>
-->
</td>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Files</strong><ul><li><a href="./Ex9.analyze.MomentCurvature2D.tcl" class="internal" title="Ex9.analyze.MomentCurvature2D.tcl">Ex9.analyze.MomentCurvature2D.tcl</a></li>
<li><a href="./MomentCurvature2D.tcl" class="internal" title="MomentCurvature2D.tcl">MomentCurvature2D.tcl</a></li></ul>
<!--
<table style="width:100%; vertical-align:top;background:#white;">
<tbody><tr><td></td></tr></tbody></table>
-->
</td>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Notes</strong>
<!--
<table style="width:100%; vertical-align:top;background:#white;">
<tbody><tr><td></td></tr></tbody></table>
-->
</td></tr></tbody></table>
<table style="margin:0; background:none; border:3px solid #ccc;">
<tbody><tr>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<h4><span class="mw-headline" id="3D_Moment-Curvature_Analysis">3D Moment-Curvature Analysis</span></h4>
<!--
<table style="width:100%; vertical-align:top;background:#white;">
  <tbody><tr><td style="color:#000;"></td></tr></tbody>
</table>
-->
</td>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Files</strong><ul><li><a href="./Ex9.analyze.MomentCurvature3D.tcl" class="internal" title="Ex9.analyze.MomentCurvature3D.tcl">Ex9.analyze.MomentCurvature3D.tcl</a></li>
<li><a href="./MomentCurvature3D.tcl" class="internal" title="MomentCurvature3D.tcl">MomentCurvature3D.tcl</a></li></ul>
<!--
<table style="width:100%; vertical-align:top;background:#white;">
  <tbody><tr><td></td></tr></tbody>
</table>
</td>
<td style="margin:0; width:25%; background:#white; vertical-align:top;">
<strong>Notes</strong>
<table style="width:100%; vertical-align:top;background:#white;">
  <tbody><tr><td></td></tr></tbody>
</table>
</td></tr>
-->
</tbody></table>
</p><p><br>
</p>


<h2><span class="mw-headline" id="Run">Run</span></h2>
<p>The model and analysis combinations for this example are numerous. The following are an small subset, for demonstration purposes:
</p>
<ul><li>To run Uniaxial-Section Model, 2D</li></ul>
<blockquote><div class="mw-highlight mw-content-ltr" dir="ltr"><pre><span></span><span class="nb">puts</span> <span class="s2">" --------------------------------- 2D Model ---------------"</span>
<span class="nb">puts</span> <span class="s2">" a. Uniaxial Section"</span>
<span class="nb">source</span> Ex9a.build.UniaxialSection2D.tcl
<span class="nb">source</span> Ex9.analyze.MomentCurvature2D.tcl
</pre></div></blockquote>
<ul><li>To run RC Section: Rectangular, Confined, Symmetric Model, 3D</li></ul>
<blockquote><div class="mw-highlight mw-content-ltr" dir="ltr"><pre><span></span><span class="nb">puts</span> <span class="s2">" --------------------------------- 3D Model ---------------"</span>
<span class="nb">puts</span> <span class="s2">" d. RC Section: Rectangular, Confined, Symmetric"</span>
<span class="nb">source</span> Ex9d.build.RCSection.RectConfinedSymm3D.tcl
<span class="nb">source</span> Ex9.analyze.MomentCurvature3D.tcl
</pre></div></blockquote>

<!--
<p>Return to <a href="/wiki/index.php/OpenSees_Examples_Manual_--_Structural_Models_%26_Analyses" class="mw-redirect" title="OpenSees Examples Manual -- Structural Models &amp; Analyses">OpenSees Examples Manual -- Structural Models &amp; Analyses</a>
</p><p>Return to <a href="/wiki/index.php/OpenSees_User" title="OpenSees User">OpenSees User</a>
</p>
-->
</div>
