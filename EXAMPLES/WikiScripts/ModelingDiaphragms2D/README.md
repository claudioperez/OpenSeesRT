# Discovering OpenSees -- Modeling Diaphragms in 2D Models with Linear and Nonlinear Elements

<p>Floor diaphragm need to be axially rigid to assure proper
distribution of seismic forces to all lateral force resisting elements
(columns and walls). Common modeling approach for frame structures is to
represent the structural components (beams and columns) by line elements
located at the original member centerlines and having cross-sectional
properties equal to those of components. The effect of a rigid diaphragm
at a floor level is usually modeled by imposing rigid constraints on all
nodes of that floor level and thus enforcing the same lateral
displacements of all nodes at the floor level. These constrains enforce
condition of zero axial strain on elements that are part of that floor
level. For sections where the neutral axis does not shift as a
consequence of bending in the beam, axial strains at element’s
centerlines are zero, and thus rigid constraints can be applied to model
the rigid diaphragm (e.g. steel sections, elastic materials). However,
for nonlinear beam-column elements (force-based or displacement-based)
with reinforced concrete fiber section where the neutral axis shifts due
to bending in the beam axial strains at element’s centerlines are no
longer zero. Thus, rigid constraints that enforce condition of zero
axial strain on elements will change the response of the frame. In this
web-learning session, the effect of a rigid constraints is demonstrated
on a 2D frame with force-based beam column elements considering three
types of sections: elastic, nonlinear steel, and nonlinear reinforced
concrete.</p>
<p>This web-learning series covers:</p>
<ul>
<li>Introduction to problem</li>
<li>Consequences of applying a rigid constraint on a force-based
beam-column element with different types of sections (demonstrated on
examples)</li>
<li>Conclusions and summary</li>
</ul>
<p><strong>PPT presentation of the seminar can be found
here:</strong></p>
<ul>
<li><a href="Media:_DiaphragmsModeling.pdf" title="wikilink">Modeling
Diaphragms in 2D Models</a></li>
</ul>
<p><strong>Video of the seminar can be found here:</strong></p>
<ul>
<li><a href="http://www.youtube.com/watch?v=_Sa5-G1lh0M">Modeling
Diaphragms in 2D Models</a></li>
</ul>
<p><strong>OpenSees files used to demonstrate the effect of rigid
constraints can be found here:</strong></p>
<ul>
<li>The main file that is to be sourced from the OpenSees interpreter:
<ul>
<li><a href="ModelingDiaphragms2D.tcl"
title="wikilink">ModelingDiaphragms2D.tcl</a></li>
</ul></li>
<li>Supporting files to be stored in the same folder with the main file:
<ul>
<li><a href="RectangularRCsection2D.tcl"
title="wikilink">RectangularRCsection2D.tcl</a> (procedure for
discretizing RC rectangular section into fibers)</li>
<li><a href="WSection.tcl" title="wikilink">WSection.tcl</a> (procedure
for discretizing W steel section into fibers)</li>
<li><a href="Media:_A10000.tcl" title="wikilink">A10000.tcl</a> (ground
motion)</li>
</ul></li>
</ul>
