# LimitStateMaterialExampleDebugged

<p>__NOTOC__</p>
<p>Application of a limit state material is demonstrated on an example
of a column modeled with force-based beam-column element. Three types of
analysis are performed on a same column:</p>
<ul>
<li>pushover analysis</li>
<li>cyclic analysis</li>
<li>time history analysis</li>
</ul>
<p><strong>Static analyses: pushover and cyclic</strong></p>
<p>Main file to source to perform pushover analysis is:</p>
<ul>
<li><a href="CenterCol_basicModel_Pushover.tcl"
title="wikilink">CenterCol_basicModel_Pushover.tcl</a></li>
</ul>
<p>Main file to source to perform cyclic analysis is:</p>
<ul>
<li><a href="CenterCol_basicModel_Cyclic.tcl"
title="wikilink">CenterCol_basicModel_Cyclic.tcl</a></li>
</ul>
<p>Supporting files to be stored in the same folder with the main
file:</p>
<ol>
<li><a href="CenterColAxialSpring.tcl"
title="wikilink">CenterColAxialSpring.tcl</a></li>
<li><a href="CenterColShearSpring.tcl"
title="wikilink">CenterColShearSpring.tcl</a></li>
<li><a href="CenterColSecFiber.tcl"
title="wikilink">CenterColSecFiber.tcl</a></li>
<li><a href="tags.tcl" title="wikilink">tags.tcl</a></li>
</ol>
<p><strong>Time history analysis</strong></p>
<p>Main file to source to perform time history analysis is:</p>
<ul>
<li><a href="CenterCol_basicModel_Dynamic.tcl"
title="wikilink">CenterCol_basicModel_Dynamic.tcl</a></li>
</ul>
<p>Supporting files to be stored in the same folder with the main
file:</p>
<ol>
<li><a href="CenterColAxialSpring.tcl"
title="wikilink">CenterColAxialSpring.tcl</a></li>
<li><a href="CenterColShearSpring.tcl"
title="wikilink">CenterColShearSpring.tcl</a></li>
<li><a href="CenterColSecFiber.tcl"
title="wikilink">CenterColSecFiber.tcl</a></li>
<li><a href="tags.tcl" title="wikilink">tags.tcl</a></li>
<li><a href="Dynamic.EQ.Uniform_LimitState.tcl"
title="wikilink">Dynamic.EQ.Uniform_LimitState.tcl</a></li>
<li><a href="LibAnalysisDynamicParameters_LimitState.tcl"
title="wikilink">LibAnalysisDynamicParameters_LimitState.tcl</a></li>
<li><a href="Media:ReadSMDFile.tcl" title="wikilink">
ReadSMDFile.tcl</a></li>
<li><a href="Media:DisplayModel2D.tcl" title="wikilink">
DisplayModel2D.tcl</a></li>
<li><a href="Media:DisplayPlane.tcl" title="wikilink">
DisplayPlane.tcl</a></li>
<li><a href="Media:_TCU047-N.tcl" title="wikilink">TCU047-N.tcl</a></li>
</ol>
<hr />
<p>Original code debugged and extended to include time history analysis
by: <span style="color:blue"> Mohammad Reza Azadi Kakavand,
College of Engineering, University of Tehran, Iran </span></p>
