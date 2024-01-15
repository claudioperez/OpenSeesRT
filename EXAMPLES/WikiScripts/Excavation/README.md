# Excavation Supported by Cantilevered Sheet Pile Wall

<p>Example prepared by: <span style="color:blue"> Christopher
McGann and Pedro Arduino, University of Washington</span></p>
<hr />
<p><a href="Examples" title="wikilink"> Return to OpenSees Examples
Page</a></p>
<p>This article describes the simulation of an excavation supported by a
cantilevered sheet pile wall using OpenSees. The model considers plane
strain conditions in two-dimensions, using quadrilateral elements with a
pressure dependent constitutive model to simulate cohesionless soil, and
beam-column elements to simulate a unit width of sheet pile wall.
Beam-solid contact elements are used to model a frictional interface
between the linear beam elements and the quadrilateral soil
elements.</p>
<p>The input file used for this example is available to view and/or
download here.</p>
<ul>
<li><a href="Media:excavation.tcl"
title="wikilink">excavation.tcl</a></li>
</ul>
<p>The pre- and post-process visualization tool <a
href="http://gid.cimne.upc.es/home">GiD</a> was used to generate this
example. Further information on how GiD can be used for geotechnical
simulations in OpenSees can be found in previous <a
href="Other_Examples" title="wikilink"> practical examples</a> posted on
this wiki. The GiD tool can be downloaded from <a
href="http://gid.cimne.upc.es/download">http://gid.cimne.upc.es/download</a>.</p>
<h2 id="model_description">Model Description</h2>
<figure>
<img src="/OpenSeesRT/contrib/static/originalMesh.png"
title="Fig. 1: Undeformed FE mesh for the excavation model."
alt="Fig. 1: Undeformed FE mesh for the excavation model." />
<figcaption aria-hidden="true">Fig. 1: Undeformed FE mesh for the
excavation model.</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/excavatedMesh.png"
title="Fig. 2: Deformed FE mesh after completion of excavation."
alt="Fig. 2: Deformed FE mesh after completion of excavation." />
<figcaption aria-hidden="true">Fig. 2: Deformed FE mesh after completion
of excavation.</figcaption>
</figure>
<p>The analysis described in this example uses the finite element mesh
shown in Fig. 1. The soil domain is 10 m tall and 10.5 m wide. The sheet
pile wall considered in this example is 0.5 m wide and has a height of
10.5 m. Beam contact elements are used to create a frictional interface
between the beam and solid elements that accounts for the full
kinematics of the beam elements and allows the constant width of the
wall to be represented in the model. To represent the 0.5 m width of the
wall, the beam contact elements (described further below) maintain a
constant distance between the wall centerline and the surrounding soil
elements of 0.25 m, resulting in the gap between the blue beam elements
and the green solid elements shown in Fig. 1.</p>
<p>The excavation is simulated by incrementally removing layers of solid
elements from one side of the sheet pile wall. After each removal, the
model is analyzed for a sufficient number of steps such that equilibrium
is reached prior to the removal of the next layer of elements. The
removal process continues until an excavation depth of 5 m is achieved.
The deformed configuration of the mesh after the completion of the
excavation analysis is shown in Fig. 2. As shown, the removal of
material from one side of the sheet pile wall allows the wall to be
pushed laterally by the remaining soil on the other side. In this
figure, the displacement magnitudes are magnified 25 times in order to
clearly display the deformation pattern.</p>
<h3 id="boundary_conditions">Boundary Conditions</h3>
<p>The nodes for the solid elements, created with two translational
degrees of freedom, are fixed against lateral translation on the left
and right hand boundaries, and fixed against vertical translation along
the lower boundary. The remaining solid element nodes are left free. The
lateral extents of the soil domain were set at ten times the wall
thickness on either side of the wall-soil interface to reduce the
effects of the boundaries on the behavior of the wall-soil system.</p>
<p>The nodes for the beam elements, created with three degrees of
freedom (two translational, one rotational), are all left free with the
exception of the node at the base of the wall, which is fixed against
vertical translation only.</p>
<h3 id="material_and_element_definitions">Material and Element
Definitions</h3>
<p>The model used in this excavation analysis contains three types of
elements, solid quadrilateral elements for the soil domain,
displacement-based beam elements for the sheet pile wall, and beam
contact elements to define the frictional wall-soil interface. Three
constitutive models (one per element type) are used to define the
material behavior of the various components in this model.</p>
<h4 id="soil_elements_and_material">Soil Elements and Material</h4>
<p>Four-node quad elements are used to model the soil using the plane
strain formulation of the <a href="Quad_Element" title="wikilink">quad
element</a>. A unit thickness is defined for these elements, indicating
that the soil domain in this plane strain model represents a 1 m thick
slice of a three-dimensional soil domain. The self-weight of the soil is
considered as a body force acting on each element with a magnitude equal
to the soil unit weight.</p>
<p>For simplicity, this example considers a homogenous cohesionless soil
domain, therefore, a single material object is used for the soil
elements. The <a
href="http://opensees.berkeley.edu/OpenSees/manuals/usermanual/1551.htm">PressureDependMultiYield02</a>
nDMaterial object is used, as this constitutive model captures the
pressure dependent strength critical to modeling a cohesionless
soil.</p>
<p>The <a href="InitialStateAnalysisWrapper"
title="wikilink">InitialStateAnalysisWrapper</a> nDMaterial enables the
use of the InitialStateAnalysis feature when defining the initial
conditions for the soil model. The initial state analysis procedure, and
the role of the <a href="InitialStateAnalysisWrapper"
title="wikilink">InitialStateAnalysisWrapper</a> material object in this
procedure, is discussed further in the following sections of this
article.</p>
<h4 id="sheet_pile_wall">Sheet Pile Wall</h4>
<p>The sheet pile wall is modeled using the <a
href="Displacement-Based_Beam-Column_Element"
title="wikilink">displacement-based beam-column elements</a> with <a
href="Elastic_Section" title="wikilink">elastic sections</a> used to
define the constitutive behavior. The beam elements are placed along the
neutral axis of the sheet pile wall, and to ensure the best behavior for
the beam contact elements, the beam elements are defined such that their
nodes are centered vertically between the nodes of the adjacent soil
elements. The beam elements extend beyond the boundaries of the soil
domain to allow the elements at the ends of the wall to be in contact
with the soil.</p>
<p>The beam-solid contact approach used in this example allows for the
use of beam-column elements to model the sheet pile wall, enabling the
results to be interpreted in the context of traditional beam theory,
i.e. it is simple to extract shear and moment diagrams from the recorded
data. Figs. 3 and 4 show the types of visualizations that can be
realized using this approach in concert with GiD, showing the evolution
of the bending moment and shear force demands in the sheet pile wall
during the excavation, respectively.</p>
<figure>
<img src="momentDiagram.gif"
title="Fig. 3: Evolution of wall bending moment demand during the excavation analysis."
alt="Fig. 3: Evolution of wall bending moment demand during the excavation analysis." />
<figcaption aria-hidden="true">Fig. 3: Evolution of wall bending moment
demand during the excavation analysis.</figcaption>
</figure>
<p><img src="shearDiagram.gif"
title="Fig. 4: Evolution of wall shear force demand during the excavation analysis."
alt="Fig. 4: Evolution of wall shear force demand during the excavation analysis." />
&lt;br style="clear: both" /&gt;</p>
<h4 id="soil_wall_interface">Soil-Wall Interface</h4>
<figure>
<img src="/OpenSeesRT/contrib/static/initialContact.png"
title="Fig. 5: Contact forces after completion of initial state analysis."
alt="Fig. 5: Contact forces after completion of initial state analysis." />
<figcaption aria-hidden="true">Fig. 5: Contact forces after completion
of initial state analysis.</figcaption>
</figure>
<figure>
<img src="/OpenSeesRT/contrib/static/finalContact.png"
title="Fig. 6: Contact forces in final excavated configuration."
alt="Fig. 6: Contact forces in final excavated configuration." />
<figcaption aria-hidden="true">Fig. 6: Contact forces in final excavated
configuration.</figcaption>
</figure>
<p>The interface between the sheet pile wall and the surrounding soil
domain is defined using two-dimensional beam-solid contact elements. The
<a href="BeamContact2D" title="wikilink">BeamContact2D</a> elements link
two master nodes (from the beam element), to a slave node (from an
adjacent soil element), and work best when the initial position of the
slave node is nearly centered vertically between the two beam nodes,
though this node need only be between the beam nodes for the contact to
be maintained. These elements use Lagrange multipliers to enforce the
contact condition, therefore an additional node is required. In the
example input file, the Lagrange multiplier nodes are numbered 1001 to
1042 and are clearly labeled with comments.</p>
<p>The constitutive behavior of the interface is defined using the <a
href="ContactMaterial2D" title="wikilink">ContactMaterial2D</a>
nDMaterial object. This constitutive model uses a regularized Coulomb
frictional law allowing for sticking and frictional slip. In this
example, a friction coefficient of 0.1 is defined, this coefficient
defines the magnitude of the frictional force that develops during the
slip condition. The second parameter, set at 1000 in this example, is
the interface stiffness which exists during the stick condition.</p>
<p>The beam-solid elements readily provide the contact forces acting on
the interface, allowing the evolution of these forces during the
excavation procedure to be tracked. Examples of the visualizations
possible for these elements are shown in Figs. 5 and 6, which show the
contact state after the completion of the initial state analysis, and
after the completion of the excavation, respectively. Fig. 12, in the
Representative Results section of this article, shows the full evolution
of the contact forces acting on the interface during the excavation
analysis.</p>
<h3 id="recorders">Recorders</h3>
<p>The large number of recorders defined for this model is necessary to
capture information from each of the layers of elements (and associated
nodes) which are removed at some point during the analysis. When the
elements and nodes for these layers are removed, the associated recorder
objects must also be removed. For each of the 10 excavation lifts in the
analysis, there is a set of recorders designated here. Each of these
sets of recorders are active during all analysis steps prior to the
removal of their associated set of elements and nodes. There is an
additional set of recorders for those nodes and elements which are never
removed, and these recorders remain active for the entire analysis.</p>
<p>The recorders are removed using their individual recorder tags.
Unlike objects such as nodes and elements, where the tag is specified in
the input file, the recorders are assigned tags based upon the order
that they appear in the input file. The tags start at 1 and end at the
total number of recorder objects. For example, the first recorder
defined for excavation lift 1 is the eleventh recorder defined in the
input file, therefore, it's associated recorder tag is 11.</p>
<h3 id="initial_state_analysis_phase">Initial State Analysis Phase</h3>
<p>For numerical analysis of geotechnical problems, it is important to
begin with the appropriate initial conditions. Use of the
InitialStateAnalysis feature allows for the initial state of stress an
strain in the soil due to the application of gravitational forces to be
achieved without any associated nodal deformation. The initiation of
this feature is designated by the command InitialStateAnalysis on After
this has been designated, the gravitational analysis is conducted. In
this example, the soil constitutive model is set to consider only linear
elastic behavior for the first 4 analysis steps, then is updated to
consider elastoplastic behavior for the remainder of the analysis. These
updates are performed using the <a
href="http://opensees.berkeley.edu/OpenSees/manuals/usermanual/1559.htm">updateMaterialStage</a>
command. The beam-contact elements are updated to ignore friction during
the gravity analysis. This is done because the application of gravity
loads in the model results in differential vertical displacements
between the beam and solid elements. If friction was considered during
this process, frictional forces be generated which do not exist in an
actual setting. The contact elements are updated to consider frictional
behavior after the conclusion of the initial state analysis phase.
Updates of these elements are performed using the <a href="setParameter"
title="wikilink">setParameter</a> command. For example, the elements are
updated to consider frictional response using, setParameter -value 1
-eleRange 1001 1042 friction Once the gravitational analysis is
complete, the conclusion of the initial state analysis feature is
designated by the command InitialStateAnalysis off When this command is
called, all displacements which had occurred during the application of
gravity are reset to zero while the <a
href="InitialStateAnalysisWrapper"
title="wikilink">InitialStateAnalysisWrapper</a> object stores the
strain which existed at the end of the gravity analysis as an initial
strain. This initial strain is then added to all subsequent strains
which are computed in the model. In the first step after switching off
the initial state analysis flag, the model is in it's original
undeformed configuration, however, the strains and stresses in the soil
elements are those that were computed during the gravitational analysis
phase. In all subsequent analysis steps, the stresses and strains which
are computed will be relative to this initial state, even though the
displacements are no longer directly associated with the computed
strains.</p>
<h3 id="excavation_analysis_phase">Excavation Analysis Phase</h3>
<p>The excavation of material is simulated in this model by removing
layers of solid elements and nodes on one side of the sheet pile wall.
After each removal, the model is allowed to reach equilibrium under the
new configuration over four analysis steps. The gravitational body
forces of the solid elements are the only loads which are applied to the
model. Removing sets of elements from one side of the wall results in an
unbalanced application of forces to the wall and the associated
deformation of the wall and surrounding soil.</p>
<p>The removal process for each excavation lift begins with the removal
of the associated recorder objects. Failure to remove the recorders
associated with removed elements and nodes will result in a segmentation
fault, so it is important to include this step. The soil elements are
removed next, followed by the contact element and Lagrange multiplier
node connected to the removed layer, and then by the nodes along the
upper boundary of the removed layer. The full removal and analysis
process for the first excavation lift is provided here</p>
<p>
```tcl
</p>
<ol>
<li>remove objects associated with lift
1-----------------------------------</li>
<li>recorders</li>
</ol>
<p>set recCount 10 for {set k 0} {$k &lt;= 8} {incr k 1} { remove
recorder [expr $recCount + $k] } set recCount [expr $recCount + $k]</p>
<ol>
<li>soil elements</li>
</ol>
<p>for {set k 1} {$k &lt;= 10} {incr k 1} { remove element [expr 190+$k]
}</p>
<ol>
<li>contact element</li>
</ol>
<p>remove element 1042</p>
<ol>
<li>lagrange multiplier node</li>
</ol>
<p>remove node 1042</p>
<ol>
<li>soil nodes</li>
</ol>
<p>remove node 430 remove node 437 remove node 446 remove node 455
remove node 461 remove node 468 remove node 473 remove node 476 remove
node 480 remove node 482 remove node 484</p>
<ol>
<li>run analysis after object removal</li>
</ol>
<p>analyze 4 
```
</p>
<p>There are a total of ten excavation lifts in the example analysis,
representing the soil on the right-hand side of the sheet pile wall up
to a depth of 5 m. Due to the use of GiD to create the input file, the
nodes are not numbered in a useful order, therefore, loops are not
employed in their removal.</p>
<h2 id="representative_results">Representative Results</h2>
<p><img src="/OpenSeesRT/contrib/static/wallMoment.png"
title="Fig. 7: Bending moment diagram for sheet pile wall after excavation. (units are kNm)"
alt="Fig. 7: Bending moment diagram for sheet pile wall after excavation. (units are kNm)" />
<img src="/OpenSeesRT/contrib/static/wallShear.png"
title="Fig. 8: Shear force diagram for sheet pile wall after excavation. (units are kN)"
alt="Fig. 8: Shear force diagram for sheet pile wall after excavation. (units are kN)" /></p>
<p>Several sets of results are provided in this section of the article
for the purposes of verifying the proper download and implementation of
the example by the user, as well as to demonstrate some of the
post-processing capabilities of GiD. The simplest means of verification
for this example is the wall bending moment and shear force diagrams
shown in Figs. 7 and 8. These diagrams are representative of the end of
the excavation analysis.</p>
<p>Figs. 9 and 10 show the distributions of vertical and shear stresses
at the end of the excavation process, respectively. These figures show
the expected concentration of stresses which occurs near the base of the
sheet pile wall, and also indicate that the level of mesh refinement in
this example is rather coarse, though it it sufficient for the purposes
of this example.</p>
<p>Fig. 11 is an animation of the full excavation analysis, showing the
nodal deformations (magnified 25 times) with the contours indicating the
magnitude of the displacements which occur. Fig. 12 shows the evolution
of the wall-soil interface contact forces as elements are removed from
the model. As shown, these force vectors do not remain horizontal,
indicating the generation of frictional forces on the interface.</p>
<p>&lt;br style="clear: both" /&gt; <img src="/OpenSeesRT/contrib/static/vertStress.png"
title="Fig. 9: Vertical stress distribution in final excavated configuration. (units are kPa)"
alt="Fig. 9: Vertical stress distribution in final excavated configuration. (units are kPa)" />
<img src="/OpenSeesRT/contrib/static/shearStress.png"
title="Fig. 10: Shear stress distribution in final excavated configuration. (units are kPa)"
alt="Fig. 10: Shear stress distribution in final excavated configuration. (units are kPa)" />
&lt;br style="clear: both" /&gt;</p>
<p><img src="nodalDisp.gif"
title="Fig. 11: Animated deformation for excavation analysis. Contours show displacement magnitude."
alt="Fig. 11: Animated deformation for excavation analysis. Contours show displacement magnitude." />
<img src="contactForces.gif"
title="Fig. 12: Evolution of wall-soil interface contact forces during excavation analysis."
alt="Fig. 12: Evolution of wall-soil interface contact forces during excavation analysis." /></p>
<p>&lt;br style="clear: both" /&gt;</p>
<p><a href="Examples" title="wikilink"> Return to OpenSees Examples
Page</a></p>
