# Pushover Analysis of 2-Story Moment Frame

<p>Example posted by: <span style="color:blue"> Laura Eads,
Stanford University</span></p>
<hr />
<p>This example demonstrates how to perform a pushover (nonlinear
static) analysis in OpenSees using a 2-story, 1-bay steel moment
resisting frame. In the first model, the nonlinear behavior is
represented using the concentrated plasticity concept with rotational
springs. In the second model, the nonlinear behavior is represented
using the distributed plasticity concept where the plastic behavior
occurs over a finite length. The rotational behavior of the plastic
regions in both models follows a bilinear hysteretic response based on
the Modified Ibarra Krawinkler Deterioration Model (Ibarra et al. 2005,
Lignos and Krawinkler 2009, 2010). For this example, all modes of cyclic
deterioration are neglected. A leaning column carrying gravity loads is
linked to the frame to simulate P-Delta effects.</p>
<p>The files needed to analyze this structure in OpenSees are included
here:</p>
<ul>
<li>The main files: <a href="pushover_concentrated.tcl"
title="wikilink">pushover_concentrated.tcl</a> and <a
href="pushover_distributed.tcl"
title="wikilink">pushover_distributed.tcl</a></li>
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
<li><a href="rotSect2DModIKModel.tcl"
title="wikilink">rotSect2DModIKModel.tcl</a> - creates a section with
bilinear rotational response that follows the Modified Ibarra Krawinkler
Deterioration Model and an elastic axial response (used in the
distributed model)</li>
<li><a href="rotLeaningCol.tcl" title="wikilink">rotLeaningCol.tcl</a> -
creates a low-stiffness rotational spring used in a leaning column</li>
</ul>
<p>All files are available in a compressed format here: <a
href="Media:pushover_example.zip"
title="wikilink">pushover_example.zip</a></p>
<p>The rest of this example describes the models and compares their
analysis results. The OpenSees models are also compared to an equivalent
model built and analyzed using the commercial program SAP2000
(http://www.csiberkeley.com/products_SAP.html).</p>
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
<h4 id="concentrated_plasticity_model_overview">Concentrated Plasticity
Model Overview</h4>
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
<figure>
<img src="/OpenSeesRT/contrib/static/Distributed_model_sketch.png"
title="Figure 2. Schematic representation of distributed plasticity OpenSees model with element number labels and [node number] labels. Note: The springs are zeroLength elements, but their sizes are greatly exaggerated in this figure for clarity."
alt="Figure 2. Schematic representation of distributed plasticity OpenSees model with element number labels and [node number] labels. Note: The springs are zeroLength elements, but their sizes are greatly exaggerated in this figure for clarity." />
<figcaption aria-hidden="true">Figure 2. Schematic representation of
distributed plasticity OpenSees model with element number labels and
[node number] labels. Note: The springs are zeroLength elements, but
their sizes are greatly exaggerated in this figure for
clarity.</figcaption>
</figure>
<h4 id="distributed_plasticity_model_overview">Distributed Plasticity
Model Overview</h4>
<p>The 2-story, 1-bay steel moment resisting frame is modeled using <a
href="Beam_With_Hinges_Element" title="wikilink">"beam with hinges"
elements</a>. Plastic hinge regions are assigned to the ends of these
elements. The plastic regions have an elastic axial response and a <a
href="Bilin_Material" title="wikilink">bilinear</a> rotational response
based on the Modified Ibarra Krawinkler Deterioration Model. An
idealized schematic of the model is presented in Figure 2.</p>
<h4 id="features_common_to_both_models">Features Common to Both
Models</h4>
<p>A leaning column with gravity loads is linked to the frame by <a
href="Truss_Element" title="wikilink">truss elements</a> to simulate
P-Delta effects.</p>
<p>To simplify this model, panel zone contributions are neglected and
plastic hinges form at the beam-column joints, i.e., centerline
dimensions are used. For an example that explicitly models the panel
zone shear distortions and includes reduced beam sections (RBS), see <a
href="Pushover_and_Dynamic_Analyses_of_2-Story_Moment_Frame_with_Panel_Zones_and_RBS"
title="wikilink">Pushover and Dynamic Analyses of 2-Story Moment Frame
with Panel Zones and RBS</a>.</p>
<p>The units of the model are kips, inches, and seconds.</p>
<h3 id="basic_geometry">Basic Geometry</h3>
<p>The basic geometry of the frame is defined by input variables for the
bay width, height of the first story, and height of a typical (i.e. not
the first) story. These values are set as WBay = 360”, HStory1 = 180”,
and HStoryTyp = 144”. The leaning column line is located one bay width
away from the frame. In addition to the nine beam-column joint nodes,
there is one additional node for each spring, which connects the spring
to the elastic element. This makes a total of 24 nodes in the structure
for the concentrated plasticity model compared to 12 nodes in the
distributed plasticity model.</p>
<h3 id="leaning_columns_and_frame_links">Leaning Columns and Frame
Links</h3>
<p>The leaning columns are modeled as <a
href="Elastic_Beam_Column_Element" title="wikilink">elastic beam-column
elements</a>. These columns have moments of inertia and areas about two
orders of magnitude larger than the frame columns in order to represent
aggregate effect of all the gravity columns (A&lt;sub&gt;leaning
column&lt;/sub&gt; = 1,000.0 in&lt;sup&gt;2&lt;/sup&gt; and
I&lt;sub&gt;leaning column&lt;/sub&gt; = 100,000.0
in&lt;sup&gt;4&lt;/sup&gt;). The columns are connected to the
beam-column joint by <a href="zeroLength_Element"
title="wikilink">zeroLength</a> rotational spring elements with very
small stiffness values so that the columns do not attract significant
moments. These springs are created using <a href="rotLeaningCol.tcl"
title="wikilink">rotLeaningCol.tcl</a>.</p>
<p><a href="Truss_Element" title="wikilink">Truss elements</a> are used
to link the frame and leaning columns and transfer the P-Delta effect.
The trusses have areas about two orders of magnitude larger than the
frame beams in order to represent aggregate effect of all the gravity
beams (A&lt;sub&gt;truss&lt;/sub&gt; = 1,000.0
in&lt;sup&gt;2&lt;/sup&gt;) and can be assumed to be axially rigid.</p>
<h3 id="rotational_springs_concentrated_plasticity_model">Rotational
Springs (Concentrated Plasticity Model)</h3>
<p>The rotational springs capture the nonlinear behavior of the frame.
As previously mentioned, the springs in the example employ a <a
href="Bilin_Material" title="wikilink">bilinear</a> hysteretic response
based on the Modified Ibarra Krawinkler Deterioration Model. Detailed
information about this model and the modes of deterioration it simulates
can be found in Ibarra et al. (2005) and Lignos and Krawinkler (2009,
2010).</p>
<p>In this example, the <a href="zeroLength_Element"
title="wikilink">zeroLength</a> spring elements connect the elastic
frame elements to the beam-column joint nodes. The springs are created
using <a href="rotSpring2DModIKModel.tcl"
title="wikilink">rotSpring2DModIKModel.tcl</a>.</p>
<h3 id="plastic_hinge_regions_distributed_plasticity_model">Plastic
Hinge Regions (Distributed Plasticity Model)</h3>
<p>In this model the plasticity is distributed over a defined length.
The axial and flexural responses of each plastic hinge region are
defined as separate <a href="Section_Command"
title="wikilink">sections</a>. The axial response is defined by an <a
href="Elastic_Section" title="wikilink">elastic section</a> while the
flexural response is defined by a <a href="Uniaxial_Section"
title="wikilink">uniaxial section</a> using the <a href="Bilin_Material"
title="wikilink">bilinear material</a> based on the Modified Ibarra
Krawinkler Deterioration Model. The responses are combined into a single
<a href="Section_Command" title="wikilink">section</a> using the <a
href="Section_Aggregator" title="wikilink">section aggregator</a>
command. P-M interactions are neglected in this example as the axial
loads are relatively low and are not expected to have a significant
influence.</p>
<h3 id="plasticity_features_common_to_both_models">Plasticity Features
Common to Both Models</h3>
<p>The input parameters for the rotational behavior of the plastic
hinges in both models are determined using empirical relationships
developed by Lignos and Krawinkler (2010) which are derived from an
extensive database of steel component tests. Alternatively, these input
parameters can be determined using approaches similar to those described
in FEMA 356 (http://www.fema.gov/library/viewRecord.do?id=1427), ATC-72
and ATC-76
(http://www.atcouncil.org/index.php?option=com_content&amp;view=article&amp;id=45&amp;Itemid=54).
In order to simplify the model, cyclic deterioration was ignored. This
was accomplished by setting all of the “L” deterioration parameter
variables to 1000.0, all of the “c” exponent variables to 1.0, and both
“D” rate of cyclic deterioration variables to 1.0.</p>
<h3
id="stiffness_modifications_to_elastic_frame_elements_concentrated_plasticity_model">Stiffness
Modifications to Elastic Frame Elements (Concentrated Plasticity
Model)</h3>
<p>Since a frame member is modeled as an elastic element connected in
series with rotational springs at either end, the stiffness of these
components must be modified so that the equivalent stiffness of this
assembly is equivalent to the stiffness of the actual frame member.
Using the approach described in Appendix B of Ibarra and Krawinkler
(2005), the rotational springs are made “n” times stiffer than the
rotational stiffness of the elastic element in order to avoid numerical
problems and allow all damping to be assigned to the elastic element. To
ensure the equivalent stiffness of the assembly is equal to the
stiffness of the actual frame member, the stiffness of the elastic
element must be “(n+1)/n” times greater than the stiffness of the actual
frame member. In this example, this is accomplished by making the
elastic element’s moment of inertia “(n+1)/n” times greater than the
actual frame member’s moment of inertia.</p>
<h3
id="stiffness_modifications_to_the_plastic_hinges_both_models">Stiffness
Modifications to the Plastic Hinges (Both Models)</h3>
<p>In order to make the nonlinear behavior of the assembly match that of
the actual frame member, the strain hardening coefficient (the ratio of
post-yield stiffness to elastic stiffness) of the plastic hinge must be
modified. If the strain hardening coefficient of the actual frame member
is denoted $\alpha$&lt;sub&gt;s,mem&lt;/sub&gt;
and the strain hardening coefficient of the spring (or plastic hinge
region) is denoted
 $\alpha_\textrm{s,spring}$  then $\alpha_\textrm{s,spring}$  $\alpha_\textrm{s,spring}$  $\alpha_\textrm{s,spring}$ 
 $\alpha_\textrm{s,spring}$  = $\alpha_\textrm{s,spring}$  $\alpha_\textrm{s,spring}$ 
&lt;math&gt;\alpha&lt;/math&gt;&lt;sub&gt;s,mem&lt;/sub&gt; / (1 + n*(1
- &lt;math&gt;\alpha&lt;/math&gt;&lt;sub&gt;s,mem&lt;/sub&gt;))</p>
<p>Note that this is a corrected version of Equation B.5 from Ibarra and
Krawinkler (2005).</p>
<h3 id="constraints">Constraints</h3>
<p>The frame columns are fixed at the base, and the leaning column is
pinned at the base. To simulate a rigid diaphragm, the horizontal
displacements of all nodes in a given floor are constrained to the
leftmost beam-column joint node using the <a href="EqualDOF_command"
title="wikilink">equalDOF</a> command.</p>
<h3 id="masses">Masses</h3>
<p>The mass is concentrated at the beam-column joints of the frame, and
each floor mass is distributed equally among the frame nodes. The mass
is assigned using the <a href="node_command" title="wikilink">node</a>
command, but it could also be assigned with the <a href="Mass_Command"
title="wikilink">mass</a> command.</p>
<h3 id="loading">Loading</h3>
<p>Gravity loads are assigned to the beam-column joint nodes using the
<a href="NodalLoad_Command" title="wikilink">nodal load command</a>.
Gravity loads tributary to the frame members are assigned to the frame
nodes while the remaining gravity loads are applied to the leaning
columns. The gravity loads are applied as a <a href="Plain_Pattern"
title="wikilink">plain load pattern</a> with a <a
href="Constant_TimeSeries" title="wikilink">constant time series</a>
since the gravity loads always act on the structure.</p>
<p>In this example, lateral loads are distributed to the frame using the
methodology of ASCE 7-10
(http://www.asce.org/Product.aspx?id=2147487569). Lateral loads are
applied to all the frame nodes in a given floor. A <a
href="Plain_Pattern" title="wikilink">plain load pattern</a> with a <a
href="Linear_TimeSeries" title="wikilink">linear time series</a> is used
for lateral load application so that loads increase with time.</p>
<h3 id="recorders">Recorders</h3>
<p>The <a href="Recorder_Command" title="wikilink">recorders</a> used in
this example include:</p>
<ul>
<li>The <a href="Drift_Recorder" title="wikilink">drift recorder</a> to
track the story and roof drift histories</li>
<li>The <a href="Node_Recorder" title="wikilink">node recorder</a> to
track the base shear reaction history</li>
<li>The <a href="Element_Recorder" title="wikilink">element recorder</a>
to track the element forces in the first story columns as well as the
moment and rotation histories of the springs in the concentrated
plasticity model</li>
</ul>
<p>To record the moment and rotation histories in the springs, the <a
href="Region_Command" title="wikilink">region command</a> was used to
assign all column springs to one group and all beam springs to a
separate group, and the region was used as an input to the <a
href="Element_Recorder" title="wikilink">element recorder</a>.</p>
<p>It is important to note that the recorders only record information
for <a href="Analyze_Command" title="wikilink">analyze commands</a> that
are called after the <a href="Recorder_Command"
title="wikilink">recorder commands</a> are called. In this example, the
recorders are placed after the gravity analysis so that the steps of the
gravity analysis do not appear in the output files.</p>
<h3 id="analysis">Analysis</h3>
<p>The structure is first analyzed under gravity loads before the
pushover analysis is conducted. The gravity loads are applied using a <a
href="Load_Control" title="wikilink">load-controlled</a> static analysis
with 10 steps. So that the gravity loads remain on the structure for all
subsequent analyses, the <a href="LoadConst_Command"
title="wikilink">loadConst command</a> is used after the gravity
analysis is completed. This command is also used to reset the time to
zero so that the pushover starts from time zero.</p>
<p>The pushover analysis is performed using a <a
href="Displacement_Control" title="wikilink">displacement-controlled</a>
static analysis. In this example, the structure was pushed to 10% roof
drift, or 32.4”. The roof node at Pier 1, node 13 in Figures 1 and 2,
was chosen as the control node where the displacement was monitored.
Incremental displacement steps of 0.01” were used. This step size was
used because it is small enough to capture the progression of hinge
formation and generate a smooth backbone curve, but not too small that
it makes the analysis time unreasonable.</p>
<h2 id="results">Results</h2>
<h3 id="comparison_of_opensees_models">Comparison of OpenSees
Models</h3>
<figure>
<img src="/OpenSeesRT/contrib/static/pushover_curve.png"
title="Figure 3. Pushover Curve: Comparison of OpenSees Models"
alt="Figure 3. Pushover Curve: Comparison of OpenSees Models" />
<figcaption aria-hidden="true">Figure 3. Pushover Curve: Comparison of
OpenSees Models</figcaption>
</figure>
<p>Theoretically, the results of the distributed plasticity model should
approach those of the concentrated plasticity model as the length of the
plastic hinge regions approaches zero. Because of localized instability
due to the stress-strain formulation of the <a
href="Beam_With_Hinges_Element" title="wikilink">beam with hinges
element</a>, the distributed plasticity model does not give reasonable
results when the length of the plastic hinge is very small (i.e.,
10e&lt;sup&gt;-5&lt;/sup&gt;). Therefore, the length of the plastic
hinges was increased from 10e&lt;sup&gt;-5&lt;/sup&gt; until the results
of this model approached those of the concentrated plasticity model. The
plastic hinge length that led to agreement between the models was 0.4%
of the frame member's total length. The periods of the concentrated and
distributed models are very close: T&lt;sub&gt;1&lt;/sub&gt; = 0.83 s
(con) vs. 0.82 s (dist) and T&lt;sub&gt;2&lt;/sub&gt; = 0.22 s (con) vs.
0.21 s (dist).</p>
<p>The results of the pushover analyses from the OpenSees models are
shown in Figure 3. This figure shows the normalized base shear (base
shear divided by the weight of the structure) versus the roof drift
(roof displacement divided by the roof elevation). The models are nearly
identical until about 2.5% roof drift when their curves begin to
diverge. The descending branch of the the concentrated plasticity model
is slightly steeper, but the two models agree reasonably well as there
is less than 10% percent difference in the base shears at 10% roof
drift.</p>
<h3 id="comparison_of_opensees_sap2000_results">Comparison of OpenSees
&amp; SAP2000 Results</h3>
<figure>
<img src="/OpenSeesRT/contrib/static/pushover_conc_plast_DRAIN_SAP.png"
title="Figure 4. Pushover Curve: Comparison OpenSees &amp; SAP2000 Models"
alt="Figure 4. Pushover Curve: Comparison OpenSees &amp; SAP2000 Models" />
<figcaption aria-hidden="true">Figure 4. Pushover Curve: Comparison
OpenSees &amp; SAP2000 Models</figcaption>
</figure>
<p>The results of the pushover analyses from the concentrated plasticity
OpenSees model and the SAP2000 model are shown in Figure 4.The OpenSees
and SAP2000 models agree very well as the difference between their base
shears at 10% roof drift is only 4%.</p>
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
<li>Lignos, D. G., and Krawinkler, H. (2012). “Sidesway Collapse of
Deteriorating Structural Systems under Seismic Excitations,” Technical
Report 177, The John A. Blume Earthquake Engineering Research Center,
Department of Civil Engineering, Stanford University, Stanford, CA.
[electronic version: <a
href="https://blume.stanford.edu/tech_reports">https://blume.stanford.edu/tech_reports</a>]</li>
<li>Lignos, D. G., and Krawinkler, H. (2011). “Deterioration Modeling of
Steel Beams and Columns in Support to Collapse Prediction of Steel
Moment Frames,” ASCE, Journal of Structural Engineering, Vol. 137 (11),
1291-1302.</li>
</ol>
