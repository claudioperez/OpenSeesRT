# Dynamic Analysis of 2-Story Moment Frame

<p>Example posted by: <span style="color:blue"> Laura Eads,
Stanford University</span>; Modified: &lt;span
style="color:red"&gt; Filipe Ribeiro, Andre Barbosa (09/03/2013)
<span style="color:blue"></p>
<hr />
<p>This example demonstrates how to perform a dynamic analysis in
OpenSees using a 2-story, 1-bay steel moment resisting frame. The
structure is subjected to the Canoga Park record from the 1994
Northridge earthquake. The nonlinear behavior is represented using the
concentrated plasticity concept with rotational springs. The rotational
behavior of the plastic regions follows a bilinear hysteretic response
based on the Modified Ibarra Krawinkler Deterioration Model (Ibarra et
al. 2005, Lignos and Krawinkler 2009, 2010). For this example, all modes
of cyclic deterioration are neglected. A leaning column carrying gravity
loads is linked to the frame to simulate P-Delta effects.</p>
<p>The files needed to analyze this structure in OpenSees are included
here:</p>
<ul>
<li>The main file: <a href="MRF_2Story_Concentrated.tcl"
title="wikilink">MRF_2Story_Concentrated.tcl</a> (last update: 10 Oct
2013)</li>
</ul>
<p>Supporting procedure files</p>
<ul>
<li><a href="DisplayModel2D.tcl" title="wikilink">DisplayModel2D.tcl</a>
- displays a 2D perspective of the model</li>
<li><a href="DisplayPlane.tcl" title="wikilink">DisplayPlane.tcl</a> -
displays a plane in the model</li>
<li><a href="rotSpring2DModIKModel.tcl"
title="wikilink">rotSpring2DModIKModel.tcl</a> - creates a bilinear
rotational spring that follows the Modified Ibarra Krawinkler
Deterioration Model (used in the concentrated model)</li>
<li><a href="rotLeaningCol.tcl" title="wikilink">rotLeaningCol.tcl</a> -
creates a low-stiffness rotational spring used in a leaning column</li>
</ul>
<p>The acceleration history for the Canoga Park record</p>
<ul>
<li><a href="NR94cnp.tcl" title="wikilink">NR94cnp.tcl</a> - contains
acceleration history in units of g</li>
</ul>
<p>All files are available in a compressed format here: <a
href="Media:dynamic_example_10Oct2013.zip"
title="wikilink">dynamic_example_10Oct2013.zip</a> (last update: 10 Oct
2013)</p>
<p>The rest of this example describes the model and shows the analysis
results.</p>
<h2 id="model_description">Model Description</h2>
<figure>
<img src="/OpenSeesRT/contrib/static/2s1b-Sketch.png"
title="Figure 1. Schematic representation of concentrated plasticity OpenSees model with element number labels and [node number] labels. Note: The springs are zeroLength elements, but their sizes are greatly exaggerated in this figure for clarity."
alt="Figure 1. Schematic representation of concentrated plasticity OpenSees model with element number labels and [node number] labels. Note: The springs are zeroLength elements, but their sizes are greatly exaggerated in this figure for clarity." />
<figcaption aria-hidden="true">Figure 1. Schematic representation of
concentrated plasticity OpenSees model with element number labels and
[node number] labels. Note: The springs are zeroLength elements, but
their sizes are greatly exaggerated in this figure for
clarity.</figcaption>
</figure>
<p>The 2-story, 1-bay steel moment resisting frame is modeled with <a
href="elastic_Beam_Column_Element" title="wikilink">elastic beam-column
elements</a> connected by <a href="zeroLength_Element"
title="wikilink">zeroLength elements</a> which serve as rotational
springs to represent the structure’s nonlinear behavior. The springs
follow a <a href="Bilin_Material" title="wikilink">bilinear</a>
hysteretic response based on the Modified Ibarra Krawinkler
Deterioration Model. A leaning column with gravity loads is linked to
the frame by <a href="Truss_Element" title="wikilink">truss elements</a>
to simulate P-Delta effects. An idealized schematic of the model is
presented in Figure 1.</p>
<p>To simplify this model, panel zone contributions are neglected,
plastic hinges form at the beam-column joints, and centerline dimensions
are used. For an example that explicitly models the panel zone shear
distortions and includes reduced beam sections (RBS), see <a
href="Pushover_and_Dynamic_Analyses_of_2-Story_Moment_Frame_with_Panel_Zones_and_RBS"
title="wikilink">Pushover and Dynamic Analyses of 2-Story Moment Frame
with Panel Zones and RBS</a>.</p>
<p>For a detailed description of this model, see <a
href="Pushover_Analysis_of_2-Story_Moment_Frame"
title="wikilink">Pushover Analysis of 2-Story Moment Frame</a>.</p>
<p>The units of the model are kips, inches, and seconds.</p>
<h2 id="damping_and_the_rayleigh_command">Damping and the Rayleigh
Command</h2>
<p>This model uses Rayleigh damping which formulates the damping matrix
as a linear combination of the mass matrix and stiffness matrix:
<strong>c</strong> = a&lt;sub&gt;0&lt;/sub&gt;*<strong>m</strong> +
a&lt;sub&gt;1&lt;/sub&gt;*<strong>k</strong>, where
a&lt;sub&gt;0&lt;/sub&gt; is the mass proportional damping coefficient
and a&lt;sub&gt;1&lt;/sub&gt; is the stiffness proportional damping
coefficient. A damping ratio of 2%, which is a typical value for steel
buildings, is assigned to the first two modes of the structure. The <a
href="Rayleigh_Damping_Command" title="wikilink">rayleigh command</a>
allows the user to specify whether the initial, current, or last
committed stiffness matrix is used in the damping matrix formulation. In
this example, only the initial stiffness matrix is used, which is
accomplished by assigning values of 0.0 to the other stiffness matrix
coefficients.</p>
<p>To properly model the structure, stiffness proportional damping is
applied only to the frame elements and not to the highly rigid truss
elements that link the frame and leaning column, nor to the leaning
column itself. OpenSees does not apply stiffness proportional damping to
<a href="zeroLength_Element" title="wikilink">zeroLength elements</a>.
In order to apply damping to only certain elements, the <a
href="Rayleigh_Damping_Command" title="wikilink">rayleigh command</a> is
used in combination with the <a href="Region_Command"
title="wikilink">region command</a>. As noted in the <a
href="Region_Command" title="wikilink">region command</a> documentation,
the region cannot be defined by BOTH elements and nodes. Because mass
proportional damping assigns damping to nodes with mass, OpenSees will
ignore any mass proportional damping that is assigned using the <a
href="Rayleigh_Damping_Command" title="wikilink">rayleigh command</a> in
combination with the <a href="Region_Command" title="wikilink">region
command</a> for a region of elements. Therefore, if using the region
command to assign damping, the mass proportional damping and stiffness
proportional damping must be assigned in separate steps.</p>
<h3
id="modifications_to_the_stiffness_proportional_damping_coefficient">Modifications
to the Stiffness Proportional Damping Coefficient</h3>
<p>As described in the “Stiffness Modifications to Elastic Frame
Elements” section of <a href="Pushover_Analysis_of_2-Story_Moment_Frame"
title="wikilink">Pushover Analysis of 2-Story Moment Frame</a>, the
stiffness of the elastic frame elements has been modified. As explained
in Ibarra and Krawinkler (2005) and Zareian and Medina (2010), the
stiffness proportional damping coefficient that is used with these
elements must also be modified. As the stiffness of the elastic elements
was made “(n+1)/n” times greater than the stiffness of the actual frame
member, the stiffness proportional damping coefficient of these elements
must also be made “(n+1)/n” times greater than the traditional stiffness
proportional damping coefficient.</p>
<h2 id="dynamic_analysis">Dynamic Analysis</h2>
<h3 id="recorders">Recorders</h3>
<p>The <a href="Recorder_Command" title="wikilink">recorders</a> used in
this example include:</p>
<ul>
<li>The <a href="Drift_Recorder" title="wikilink">drift recorder</a> to
track the story and roof drift histories</li>
<li>The <a href="Node_Recorder" title="wikilink">node recorder</a> to
track the floor displacement and base shear reaction histories</li>
<li>The <a href="Element_Recorder" title="wikilink">element recorder</a>
to track the element forces in the first story columns as well as the
moment and rotation histories of the springs in the concentrated
plasticity model</li>
</ul>
<p>For the <a href="Element_Recorder" title="wikilink">element
recorder</a>, the <a href="Region_Command" title="wikilink">region
command</a> was used to assign all column springs to one group and all
beam springs to a separate group.</p>
<p>It is important to note that the recorders only record information
for <a href="Analyze_Command" title="wikilink">analyze commands</a> that
are called after the <a href="Recorder_Command"
title="wikilink">recorder commands</a> are called. In this example, the
recorders are placed after the gravity analysis so that the steps of the
gravity analysis do not appear in the output files.</p>
<h3 id="analysis">Analysis</h3>
<p>The structure is analyzed under gravity loads before the dynamic
analysis is conducted. The gravity loads are applied using a <a
href="Load_Control" title="wikilink">load-controlled</a> static analysis
with 10 steps. So that the gravity loads remain on the structure for all
subsequent analyses, the <a href="LoadConst_Command"
title="wikilink">loadConst command</a> is used after the gravity
analysis is completed. This command is also used to reset the time to
zero so that the dynamic analysis starts from time zero.</p>
<p>For the dynamic analysis, the structure is subjected to the Canoga
Park record from the 1994 Northridge earthquake. To apply the ground
motion to the structure, the <a href="Uniform_Exciatation_Pattern"
title="wikilink">uniform excitation pattern</a> is used. The name of the
file containing the acceleration record, timestep of the ground motion,
scale factor applied to the ground motion, and the direction in which
the motion is to be applied must all be specified as part of the <a
href="Uniform_Exciatation_Pattern" title="wikilink">uniform excitation
pattern command</a>.</p>
<p>To execute the dynamic analysis, the <a href="Analyze_Command"
title="wikilink">analyze command</a> is used with the specified number
of analysis steps and the timestep of the analysis. The timestep used in
the analysis should be less than or equal to the timestep of the input
ground motion.</p>
<h2 id="results">Results</h2>
<figure>
<img src="/OpenSeesRT/contrib/static/Dhist_plot_ConcDynam.png"
title="Figure 2. Floor Displacement History"
alt="Figure 2. Floor Displacement History" />
<figcaption aria-hidden="true">Figure 2. Floor Displacement
History</figcaption>
</figure>
<p>The floor displacement histories from the dynamic analysis are shown
in Figure 2. The top graph shows the ground acceleration history while
the middle and bottom graphs show the displacement time histories of the
3rd floor (roof) and 2nd floor, respectively.</p>
<h2 id="references">References</h2>
<ol>
<li>Ibarra, L. F., and Krawinkler, H. (2005). “Global collapse of frame
structures under seismic excitations,” Technical Report 152, The John A.
Blume Earthquake Engineering Research Center, Department of Civil
Engineering, Stanford University, Stanford, CA. [electronic version: <a
href="https://blume.stanford.edu/tech_reports">https://blume.stanford.edu/tech_reports</a>]</li>
<li>Ibarra, L. F., Medina, R. A., and Krawinkler, H. (2005). “Hysteretic
models that incorporate strength and stiffness deterioration,”
Earthquake Engineering and Structural Dynamics, Vol. 34, 12, pp.
1489-1511.</li>
<li>Lignos, D. G., and Krawinkler, H. (2009). “Sidesway Collapse of
Deteriorating Structural Systems under Seismic Excitations,” Technical
Report 172, The John A. Blume Earthquake Engineering Research Center,
Department of Civil Engineering, Stanford University, Stanford, CA.</li>
<li>Lignos, D. G., and Krawinkler, H. (2011). “Deterioration Modeling of
Steel Beams and Columns in Support to Collapse Prediction of Steel
Moment Frames,” ASCE, Journal of Structural Engineering, Vol. 137 (11),
1291-1302.</li>
<li>Zareian, F. and Medina, R. A. (2010). “A practical method for proper
modeling of structural damping in inelastic plane structural systems,”
Computers &amp; Structures, Vol. 88, 1-2, pp. 45-53.</li>
</ol>
