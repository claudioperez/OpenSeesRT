# Dynamic Analyses of 1-Story Moment Frame with Viscous Dampers

Example posted by:
***[Sarven Akcelyan](http://sarvenakcelyan.com)***
 and 
***[Prof. Dimitrios G. Lignos](http://dimitrios-lignos.research.mcgill.ca/PLignos.html)***
(McGill University)

------------------------------------------------------------------------

This example demonstrates how to use the viscous damper material within
a simple single story shear frame.

The files needed to analyze this structure in OpenSees are included
here:

-   The main file:

Supporting files

-   TakY.th - uses the JR Takatori record from the Kobe 1995 earthquake
    (available in the zip file below)

All files are available in a compressed format here:
[Viscous_Damper_Example.zip](Media:Viscous_Damper_Example.zip "wikilink")

The rest of this example describes the model and shows the analysis
results.

## Model Description

<figure>
<img src="Viscous-Fig1.jpg"
title="Figure 1. Schematic representation of a viscous damper installed in the single story moment resisting frame."
width="300"
alt="Figure 1. Schematic representation of a viscous damper installed in the single story moment resisting frame." />
<figcaption aria-hidden="true">Figure 1. Schematic representation of a
viscous damper installed in the single story moment resisting
frame.</figcaption>
</figure>

The viscous damper is modeled with a [Two Node Link
Element](Two_Node_Link_Element "wikilink"). This element follows a
[Viscous damper](ViscousDamper_Material "wikilink") hysteretic response.
An idealized schematic of the model is presented in Figure 1.

The units of the model are mm, kN, and seconds.

### Basic Geometry 

The single bay single story frame shown in Figure 1 has 5000mm bay width
and 3000mm story height (centerline). The period of the system is
**0.7sec**. Columns and beams of the frame are modeled with elastic
beam-column elements.

### Damper Links

[A Two Node Link Element](Two_Node_Link_Element "wikilink") is used to
link the two nodes that define the geometry of the viscous damper.

### Constraints

the Nodes at the base of the frame are fixed. The beam (element 3 in
Figure 1) is considered to be rigid.

### Viscous Damper Material

To model the viscous damper the
[ViscousDamper](ViscousDamper_Material "wikilink") is used. The input
parameters that are selected for the damper example are as follows:
Axial Stiffness K = 25 kN/mm, Damping Coefficient Cd=20.74
kN(s/mm)\<sup\>0.35\</sup\> and exponent a=0.35.

### Loading

The single story frame with viscous damper is subjected to the 50% JR
Takatori record from the Kobe 1995 earthquake in Japan.

### Recorders

The [recorders](Recorder_Command "wikilink") used in this example
include:

-   The [Element recorder](Element_Recorder "wikilink") to track the
    damper axial force and axial displacement.
-   The [Node recorder](Node_Recorder "wikilink") to track the Frame
    displacement history at its roof.

### Analysis

A uniform excitation option is selected with application of ground
acceleration history as the imposed motion. The Newmark integration
scheme is selected for integration of the equations of motion with a
time step dt = 0.001sec. Two percent mass proportional damping is used.

## Results

### Simulation Results for the 50% JR Takatori Record

<figure>
<img src="./Viscous-Fig2.png"
title="Figure 2. Displacement history at the roof of the single story MRF"
width="700"
alt="Figure 2. Displacement history at the roof of the single story MRF" />
<figcaption aria-hidden="true">Figure 2. Displacement history at the
roof of the single story MRF</figcaption>
</figure>

-   The force - displacement relationship from the viscous damper are
    shown in Figure 3. A comparison with results from a SAP2000 model is
    also shown in Figure 3. Results are nearly identical between the two
    models.

<figure>
<img src="./Viscous1-Fig3.png"
title="Figure 3. Force - displacement relationship of the viscous damper and comparison with identical model in SAP2000"
width="500"
alt="Figure 3. Force - displacement relationship of the viscous damper and comparison with identical model in SAP2000" />
<figcaption aria-hidden="true">Figure 3. Force - displacement
relationship of the viscous damper and comparison with identical model
in SAP2000</figcaption>
</figure>
