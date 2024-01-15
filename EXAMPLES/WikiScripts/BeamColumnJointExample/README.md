# [BeamColumnJointExample Corrected](https://opensees.berkeley.edu/wiki/index.php/BeamColumnJointExample_Corrected)

<p>The example files (<a href="Media:PR.tcl"
title="wikilink">PR.tcl</a>, <a href="Media:procMKPC2.tcl"
title="wikilink">procMKPC2.tcl</a>, <a
href="Media:procUniaxialPinching.tcl"
title="wikilink">procUniaxialPinching.tcl</a>, <a
href="Media:procRC2.tcl" title="wikilink">procRC2.tcl</a>) create a
model of a RC beam column sub-assemblage. The cruciform is subjected to
constant gravity load at nodes 4 and 7 and pseudo-static cyclic lateral
load under displacement control at node 10. The beam-column-joint region
(element number 7) is represented using a beamColumnJoint element, and
the beams and columns (element numbers 1 through 6) are modeled using
the nonlinearBeamColumn element. The beam-column joint consists of 13
components that may have different material constitutive models; in this
example 9 of the 13 components utilize the nonlinear material model -
Pinching4.</p>
<figure>
<img src="BeamColumnJointExample.jpg" title="BeamColumnJointExample.jpg"
alt="BeamColumnJointExample.jpg" />
<figcaption aria-hidden="true">BeamColumnJointExample.jpg</figcaption>
</figure>
<p>The displacement history for node 10 is as shown below.</p>
<figure>
<img src="BeamColumnJointExample2.jpg"
title="BeamColumnJointExample2.jpg" alt="BeamColumnJointExample2.jpg" />
<figcaption aria-hidden="true">BeamColumnJointExample2.jpg</figcaption>
</figure>
<p>The p-delta response of cruciform, along with the response of each of
the nonlinear joint components is shown below.</p>
<figure>
<img src="BeamColumnJointExample3.jpg"
title="BeamColumnJointExample3.jpg" alt="BeamColumnJointExample3.jpg" />
<figcaption aria-hidden="true">BeamColumnJointExample3.jpg</figcaption>
</figure>
<p>The shear panel response shows the moment-curvature relationships
whereas the bar slips at the beam top and bottom are represented by the
force-slip plots</p>
<figure>
<img src="BeamColumnJointExample4.jpg"
title="BeamColumnJointExample4.jpg" alt="BeamColumnJointExample4.jpg" />
<figcaption aria-hidden="true">BeamColumnJointExample4.jpg</figcaption>
</figure>
