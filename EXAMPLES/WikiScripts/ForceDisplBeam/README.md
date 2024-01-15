# Discovering OpenSees -- Force-based Element vs. Displacement-based Element

<p>The two most commonly used OpenSees elements for modeling beam-column
elements are force-based element (FBE) and displacement-based element
(DBE). Although OpenSees command for defining these two elements has the
same arguments, a beam-column element needs to be modeled differently
using these two elements to achieve a comparable level of accuracy. The
intent of this seminar is to show users how to properly model
beam-column elements with both FBE and DBE. The theory behind these two
elements along with two examples will be presented to users to enhance
their understanding of the elements and assure their correct
application. This seminar covers:</p>
<ul>
<li>Introduction</li>
<li>Theory of force-based element</li>
<li>Theory of displacement-based element</li>
<li>Example 1 - steel beam (made up configuration that highlights
modeling differences between FBE and DBE elements)</li>
<li>Example 2 - practical example (bridge RC column)</li>
<li>Summary and conclusions</li>
</ul>
<p><strong>PPT presentation of the seminar can be found
here:</strong></p>
<ul>
<li><a href="Media:_FBEvsDBE_final.pdf" title="wikilink">FBE vs.
DBE</a></li>
</ul>
<p><strong>Video of the seminar can be found here:</strong></p>
<ul>
<li><a href="http://www.youtube.com/watch?v=yk-1k2aF53E">FBE vs.
DBE</a></li>
</ul>
<p><strong>OpenSees files used to demonstrate the effect of rigid
constraints can be found here:</strong></p>
<ul>
<li>The main file that is to be sourced from the OpenSees interpreter:
<ul>
<li><a href="ConventionalColumn_Cyclic.tcl"
title="wikilink">ConventionalColumn_Cyclic.tcl</a></li>
</ul></li>
<li>Supporting files to be stored in the same folder with the main file:
<ul>
<li><a href="LibUnits.tcl" title="wikilink">LibUnits.tcl</a> (define
system of units)</li>
<li><a href="singlecycle.tcl" title="wikilink">singlecycle.tcl</a>
(procedure for writing one cycle of displacement history)</li>
<li><a href="Media:_leh415.xls" title="wikilink">leh415.xls</a>
(experimental force-displacement response)</li>
</ul></li>
</ul>
