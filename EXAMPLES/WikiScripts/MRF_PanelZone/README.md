# Pushover and Dynamic Analyses of 2-Story Moment Frame with Panel Zones and RBS

<p>Example posted by: <span style="color:blue"> Laura Eads,
Stanford University</span></p>
<hr />
<p>This example is an extension of the <a
href="Pushover_Analysis_of_2-Story_Moment_Frame"
title="wikilink">Pushover Analysis of 2-Story Moment Frame</a> and <a
href="Dynamic_Analysis_of_2-Story_Moment_Frame" title="wikilink">
Dynamic Analysis of 2-Story Moment Frame</a> examples which illustrates
the explicit modeling of shear distortions in panel zones and uses
reduced beam sections (RBS) which are offset from the panel zones. Both
pushover and dynamic analyses are performed in this example. The
structure is the same 2-story, 1-bay steel moment resisting frame used
in the other examples where the nonlinear behavior is represented using
the concentrated plasticity concept with rotational springs. The
rotational behavior of the plastic regions follows a bilinear hysteretic
response based on the Modified Ibarra Krawinkler Deterioration Model
(Ibarra et al. 2005, Lignos and Krawinkler 2009, 2010). For this
example, all modes of cyclic deterioration are neglected. A leaning
column carrying gravity loads is linked to the frame to simulate P-Delta
effects.</p>
<p>The files needed to analyze this structure in OpenSees are included
here:</p>
<ul>
<li>The main file: <a href="MRF_2Story_Concentrated_PanelZone.tcl"
title="wikilink">MRF_2Story_Concentrated_PanelZone.tcl</a></li>
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
<li><a href="rotPanelZone2D.tcl" title="wikilink">rotPanelZone2D.tcl</a>
- creates a rotational spring to capture panel zone shear
distortions</li>
<li><a href="elemPanelZone2D.tcl"
title="wikilink">elemPanelZone2D.tcl</a> -creates eight elastic elements
which form a rectangular panel zone</li>
</ul>
<p>The acceleration history for the Canoga Park record</p>
<ul>
<li><a href="NR94cnp.tcl" title="wikilink">NR94cnp.tcl</a> - contains
acceleration history in units of g</li>
</ul>
<p>All files are available in a compressed format here: <a
href="Media:MRF_PanelZone_example.zip" title="wikilink">
MRF_PanelZone_example.zip</a></p>
<p>The rest of this example describes the model and presents the
analysis results. The OpenSees model is also compared to an equivalent
model built and analyzed using the commercial program SAP2000
(http://www.csiberkeley.com/products_SAP.html).</p>
<h2 id="model_description">Model Description</h2>
<figure>
<img src="/OpenSeesRT/contrib/static/_Eads_conc_PanelZone_sketch.png"
title="Figure 1. Schematic representation of concentrated plasticity OpenSees model including explicit modeling of the panel zones. Element number labels and [node number] labels are also shown. A detailed view of a typical panel zone is presented in Figure 2. Note: The springs are zeroLength elements, but their sizes are greatly exaggerated in this figure for clarity."
alt="Figure 1. Schematic representation of concentrated plasticity OpenSees model including explicit modeling of the panel zones. Element number labels and [node number] labels are also shown. A detailed view of a typical panel zone is presented in Figure 2. Note: The springs are zeroLength elements, but their sizes are greatly exaggerated in this figure for clarity." />
<figcaption aria-hidden="true">Figure 1. Schematic representation of
concentrated plasticity OpenSees model including explicit modeling of
the panel zones. Element number labels and [node number] labels are also
shown. A detailed view of a typical panel zone is presented in Figure 2.
Note: The springs are zeroLength elements, but their sizes are greatly
exaggerated in this figure for clarity.</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/_Eads_PanelZone_NumberingConventions.png"
title="Figure 2. Schematic representation of a typical panel zone with element number labels and [node number] labels shown. Note: The spring is a zeroLength element, but its size is greatly exaggerated in this figure for clarity."
alt="Figure 2. Schematic representation of a typical panel zone with element number labels and [node number] labels shown. Note: The spring is a zeroLength element, but its size is greatly exaggerated in this figure for clarity." />
<figcaption aria-hidden="true">Figure 2. Schematic representation of a
typical panel zone with element number labels and [node number] labels
shown. Note: The spring is a zeroLength element, but its size is greatly
exaggerated in this figure for clarity.</figcaption>
</figure>
<p>The 2-story, 1-bay steel moment resisting frame is modeled with <a
href="elastic_Beam_Column_Element" title="wikilink">elastic beam-column
elements</a> connected by <a href="zeroLength_Element"
title="wikilink">zeroLength elements</a> which serve as rotational
springs to represent the structure’s nonlinear behavior. The springs
follow a <a href="Bilin_Material" title="wikilink">bilinear</a>
hysteretic response based on the Modified Ibarra Krawinkler
Deterioration Model. The panel zones are explicitly modeled with eight
<a href="elastic_Beam_Column_Element" title="wikilink">elastic
beam-column elements</a> and one <a href="zeroLength_Element"
title="wikilink">zeroLength element</a> which serves as rotational
spring to represent shear distortions in the panel zone. A leaning
column with gravity loads is linked to the frame by <a
href="Truss_Element" title="wikilink">truss elements</a> to simulate
P-Delta effects. An idealized schematic of the model is presented in
Figure 1.</p>
<p>A detailed description of this model is provided in <a
href="Pushover_Analysis_of_2-Story_Moment_Frame"
title="wikilink">Pushover Analysis of 2-Story Moment Frame</a>. This
section merely highlights the important differences in this model,
namely the inclusion of panel zones and reduced beam sections (RBS)
which are offset from the panel zones.</p>
<p>The units of the model are kips, inches, and seconds.</p>
<h3 id="panel_zones">Panel Zones</h3>
<p>The panel zone is the joint region where beams and columns intersect.
In this model it consists of the rectangular area of the column web that
lies between the flanges of the connecting beam(s). The panel zone
deforms primarily in shear due to the opposing moments in the beams and
columns. To capture these deformations, the panel zone is explicitly
modeled using the approach of Gupta and Krawinkler (1999) as a rectangle
composed of eight very stiff <a href="elastic_Beam_Column_Element"
title="wikilink">elastic beam-column elements</a> with one <a
href="zeroLength_Element" title="wikilink">zeroLength element</a> which
serves as rotational spring to represent shear distortions in the panel
zone (see Figure 2). At the three corners of the panel zone without a
spring, the elements are joined by a simple pin connection which is
achieved by using the <a href="EqualDOF_command"
title="wikilink">equalDOF</a> command to constrain both translational
degrees of freedom. The eight <a href="elastic_Beam_Column_Element"
title="wikilink">elastic beam-column elements</a> each have an area of
1,000.0 in&lt;sup&gt;2&lt;/sup&gt; and a moment of inertia equal to
10,000.0 in&lt;sup&gt;4&lt;/sup&gt; in order to give them high axial and
flexural stiffness, respectively. The elements are defined in <a
href="elemPanelZone2D.tcl" title="wikilink">elemPanelZone2D.tcl</a>. The
spring has a trilinear backbone which is created with the <a
href="Hysteretic_Material" title="wikilink">Hysteretic material</a> in
<a href="rotPanelZone2D.tcl" title="wikilink">rotPanelZone2D.tcl</a>.
This procedure also constrains the translational degrees of freedom at
the corners of the panel zone. The spring’s backbone curve is derived
using the principle of virtual work applied to a deformed configuration
of the panel zone (Gupta and Krawinkler 1999).</p>
<h3 id="reduced_beam_sections_rbs">Reduced Beam Sections (RBS)</h3>
<p>Using an RBS which is offset from the beam-column joint ensures that
the beam’s plastic hinge forms away from the column and thus protects
the column’s integrity. In this model, the decrease in moment of inertia
at the RBS is neglected; however, the yield moment at the RBS is
calculated based on the reduced section properties. The plastic hinge is
modeled by a rotational spring placed at the center of the RBS. An <a
href="elastic_Beam_Column_Element" title="wikilink">elastic beam-column
element</a> is used to connect the spring and the panel zone. Since this
element is not part of the spring-elastic element-beam subassembly
described in the “Stiffness Modifications to Elastic Frame Elements”
section of the <a href="Pushover_Analysis_of_2-Story_Moment_Frame"
title="wikilink">Pushover Analysis of 2-Story Moment Frame</a> example,
its moment of inertia and stiffness proportional damping coefficient are
not modified by an “n” factor.</p>
<h3 id="application_points_for_masses_and_loading">Application Points
for Masses and Loading</h3>
<p>Since loads cannot be applied at the center of the beam-column joint,
gravity loads are applied at the top node of the panel zone where it
meets the column (node xy7 in Figure 2). Both masses and lateral loads
are applied at the centerline of the floor level along the right side of
the panel zone (node xy05 in Figure 2).</p>
<h2 id="analysis">Analysis</h2>
<figure>
<img src="/OpenSeesRT/contrib/static/_Eads_Pushover_PZ_SAP.png"
title="Figure 3. Pushover Curve: Comparison OpenSees &amp; SAP2000 Models"
alt="Figure 3. Pushover Curve: Comparison OpenSees &amp; SAP2000 Models" />
<figcaption aria-hidden="true">Figure 3. Pushover Curve: Comparison
OpenSees &amp; SAP2000 Models</figcaption>
</figure>
<h3 id="pushover">Pushover</h3>
<p>The pushover analysis is identical to the analysis performed in the
<a href="Pushover_Analysis_of_2-Story_Moment_Frame"
title="wikilink">Pushover Analysis of 2-Story Moment Frame</a> example
where the structure is pushed to 10% roof drift, or 32.4”.</p>
<h3 id="dynamic">Dynamic</h3>
<p>The dynamic analysis is identical to the analysis performed in the <a
href="Dynamic_Analysis_of_2-Story_Moment_Frame" title="wikilink">
Dynamic Analysis of 2-Story Moment Frame</a> example where the structure
is subjected to the 1994 Northridge Canoga Park record.</p>
<h2 id="results">Results</h2>
<p>The first and second mode periods of the structure obtained from an
eigenvalue analysis are T&lt;sub&gt;1&lt;/sub&gt; = 0.81 s and
T&lt;sub&gt;2&lt;/sub&gt; = 0.18 s, respectively. These values agree
with the SAP2000 model which had periods of T&lt;sub&gt;1&lt;/sub&gt; =
0.81 s and T&lt;sub&gt;2&lt;/sub&gt; = 0.20 s.</p>
<p>The periods of this OpenSees model are slightly smaller than the
periods of the structure used in <a
href="Pushover_Analysis_of_2-Story_Moment_Frame"
title="wikilink">Pushover Analysis of 2-Story Moment Frame</a> which had
periods of T&lt;sub&gt;1&lt;/sub&gt; = 0.83 s and
T&lt;sub&gt;2&lt;/sub&gt; = 0.22 s. This is expected because the
including the panel zone regions makes the structure stiffer.</p>
<h3 id="pushover_results">Pushover Results</h3>
<figure>
<img src="/OpenSeesRT/contrib/static/_Eads_PZ_thist.png"
title="Figure 4. Acceleration and Floor Displacement Histories"
alt="Figure 4. Acceleration and Floor Displacement Histories" />
<figcaption aria-hidden="true">Figure 4. Acceleration and Floor
Displacement Histories</figcaption>
</figure>
<p>A comparison of the pushover results from the OpenSees and SAP2000
models is shown in Figure 3. As demonstrated by this figure, the results
are nearly identical.</p>
<h3 id="dynamic_results">Dynamic Results</h3>
<p>The floor displacement histories from the dynamic analysis are shown
in Figure 4. The top graph shows the ground acceleration history while
the middle and bottom graphs show the displacement time histories of the
3rd floor (roof) and 2nd floor, respectively</p>
<h2 id="references">References</h2>
<ol>
<li>Gupta, A., and Krawinkler, H. (1999). "Seismic Demands for
Performance Evaluation of Steel Moment Resisting Frame Structures,"
Technical Report 132, The John A. Blume Earthquake Engineering Research
Center, Department of Civil Engineering, Stanford University, Stanford,
CA. [electronic version: <a
href="https://blume.stanford.edu/tech_reports">https://blume.stanford.edu/tech_reports</a>]</li>
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
Steel Components in Support of Collapse Prediction of Steel Moment
Frames under Earthquake Loading", ASCE, Journal of Structural
Engineering, Vol. 137 (11), 1291-1302.</li>
</ol>
