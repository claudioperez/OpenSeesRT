# Simple Truss

The first example is a simple truss structure. The purpose of this
example is to show that model generation in OpenSees can resemble
typical finite element analysis programs with the definition of nodes,
materials, elements, loads and constraints. The example also
demonstrates how an analysis object is 'built' from component objects.

## Example 1.1 -- Linear elasticity

This example is of a linear-elastic three bar truss, as shown in
figure [example1](#example1){reference-type="ref"
reference="example1"}, subject to static loads.

1.  `Example1.1.tcl`

The model consists of four nodes, three truss elements, a single load
pattern with a nodal load acting at node 4, and constraints at the three
support nodes. Since the truss elements have the same elastic material,
a single Elastic material object is created.

::: {.center}
![Example 1.1](../fig_files/Example1.svg)
:::

The model is linear, so we use a solution Algorithm of type `Linear`. Even
though the solution is linear, we have to select a procedure for
applying the load, which is called an `Integrator`. For this problem, a
`LoadControl` integrator advances the solution. The equations are formed
using a banded system, so the System is `BandSPD` (banded, symmetric
positive definite). This is a good choice for most moderate size models.
The equations have to be numbered, so typically an RCM numberer object
is used (for Reverse Cuthill-McKee). The constraints are most easily
represented with a `Plain` constraint handler.

Once all the components of an analysis are defined, the Analysis object
itself is created. For this problem a Static Analysis object is used.

When the analysis is complete the state of node 4 and all three elements
will be printed to the screen. Nothing is recorded for later use.

The Tcl script for the example is shown below. A comment is indicated by
a `#` symbol. In the comments below, the syntax for important commands
are given.

<details class="note">
<summary>Tcl source</summary>

    # OpenSees Example 1.1
    # OpenSees Primer
    #
    # Units: kips, in, sec

    # ------------------------------
    # Start of model generation
    # ------------------------------

    # Create ModelBuilder (with two-dimensions and 2 DOF/node)
    model BasicBuilder -ndm 2 -ndf 2

    # Create nodes
    # ------------

    # Create nodes & add to Domain - command: node nodeId xCrd yCrd
    node 1   0.0  0.0
    node 2 144.0  0.0
    node 3 168.0  0.0
    node 4  72.0 96.0

    # Set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?
    fix 1 1 1 
    fix 2 1 1
    fix 3 1 1


    # Define materials for truss elements
    # -----------------------------------

    # Create Elastic material prototype - command: uniaxialMaterial Elastic matID E
    uniaxialMaterial Elastic 1 3000

    # Define elements
    # ---------------


    # Create truss elements - command: element truss trussID node1 node2 A matID
    element truss 1 1 4 10.0 1
    element truss 2 2 4 5.0 1
    element truss 3 3 4 5.0 1


    # Define loads
    # ------------

    # Create a Plain load pattern with a linear TimeSeries
    pattern Plain 1 "Linear" {
        # Create the nodal load - command: load nodeID xForce yForce
        load 4 100 -50
    }

    # End of model generation
    # ------------------------------


    # ------------------------------
    # Start of analysis generation
    # ------------------------------

    # Create the system of equation, a SPD using a band storage scheme
    system BandSPD

    # Create the DOF numberer, the reverse Cuthill-McKee algorithm
    numberer RCM

    # Create the constraint handler, a Plain handler is used as homo constraints
    constraints Plain

    # Create the integration scheme, the LoadControl scheme using steps of 1.0
    integrator LoadControl 1.0

    # Create the solution algorithm, a Linear algorithm is created
    algorithm Linear

    # create the analysis object 
    analysis Static 


    # End of analysis generation
    # ------------------------------


    # ------------------------------
    # Start of recorder generation
    # ------------------------------

    # create a Recorder object for the nodal displacements at node 4
    recorder Node -file example.out -load -node 4 -dof 1 2 disp

    # --------------------------------
    # End of recorder generation
    # ---------------------------------


    # ------------------------------
    # Finally perform the analysis
    # ------------------------------

    # Perform the analysis
    analyze 1

    # Print the current state at node 4 and at all elements
    print node 4
    print ele


</details>


     Node: 4
            Coordinates  : 72 96 
            commitDisps: 0.530093 -0.177894 
            unbalanced Load: 100 -50 

    Element: 1 type: Truss  iNode: 1 jNode: 4 Area: 10 Total Mass: 0 
             strain: 0.00146451 axial load: 43.9352 
             unbalanced load: -26.3611 -35.1482 26.3611 35.1482 
             Material: Elastic tag: 1
             E: 3000 eta: 0

    Element: 2 type: Truss  iNode: 2 jNode: 4 Area: 5 Total Mass: 0 
             strain: -0.00383642 axial load: -57.5463 
             unbalanced load: -34.5278 46.0371 34.5278 -46.0371 
             Material: Elastic tag: 1
             E: 3000 eta: 0

    Element: 3 type: Truss  iNode: 3 jNode: 4 Area: 5 Total Mass: 0 
             strain: -0.00368743 axial load: -55.3114 
             unbalanced load: -39.1111 39.1111 39.1111 -39.1111 
             Material: Elastic tag: 1
             E: 3000 eta: 0

For the node, displacements and loads are given. For the truss elements,
the axial strain and force are provided along with the resisting forces
in the global coordinate system.

The file **example.out**, specified in the recorder command, provides
the nodal displacements for the x and y directions of node 4. The file
consists of a single line of code:

    1.0 0.530093 -0.177894 

The $1.0$ corresponds to the load factor (pseudo time) in the model at
which point the recorder was invoked. The $0.530093$ and $-0.177894$
correspond to the response at node $4$ for the 1 and 2
degree-of-freedom. Note that if more analysis steps had been peformed,
the line would contain a line for every analysis step that completed
successfully.

