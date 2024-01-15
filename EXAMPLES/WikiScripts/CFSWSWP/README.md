# Cold-Formed Steel Wood Sheathed Shear Wall Panel examples

<p>In order to check the accuracy of the proposed models, results from
quasi-static non-linear analyses of Wood-sheathed Cold-Formed Steel
Shear Wall Panel have been carried out using CFSWSWP uniaxialMaterial
model implemented into OpenSees version 2.4.6. For this purpose,
specimen No 26, No 32 and No 16 tested by Branston et al. (2006) were
selected from the literature and analyzed under similar loading
conditions. This set of specimens covers a wide range of variation in
physical and mechanical characteristics such as: spacing, number, shear
strength, diameter of screw fasteners (sc, nc, Vs, and ds); wall aspect
ratio (H/W); frame thickness, chord stud moment of inertia, chord stud
cross section area, interior stud moment of inertia, yield and tensile
strengths of steel frame (tf, Af, Ife, Ifi, fyf, and fuf); type of
sheathing wood (DFP: Douglas Fir Plywood, OSB: Oriented strand Board and
CSP: Canadian Soft Plywood), sheathing thickness, as well as the yield
and tensile strengths of sheathing (type, ts, Fys, and Fus).</p>
<p>Each SWP was subjected to the Consortium of Universities for Research
in Earthquake Engineering (CUREE) loading protocol.</p>
<p><strong>TCL FILES:</strong></p>
<p><a href="Media:26n.tcl" title="wikilink">26n.tcl</a></p>
<p><a href="Media:32n.tcl" title="wikilink">32n.tcl</a></p>
<p><a href="Media:16n.tcl" title="wikilink">16n.tcl</a></p>
<p><a href="Media:procRCycDAns.tcl"
title="wikilink">procRCycDAns.tcl</a></p>
<p>uniaxialMaterial CFSWSWP 26 2440 2440 344 1.12 181600 51240 12.5 1
4.064 3256 102 98 3 0 0</p>
<figure>
<img src="/OpenSeesRT/contrib/static/26.png"
title="Fig 1. Comparison between wood-sheathed CFS SWP experimental and numerical results specimens No 26"
width="400"
alt="Fig 1. Comparison between wood-sheathed CFS SWP experimental and numerical results specimens No 26" />
<figcaption aria-hidden="true">Fig 1. Comparison between wood-sheathed
CFS SWP experimental and numerical results specimens No 26</figcaption>
</figure>
<p>uniaxialMaterial CFSWSWP 32 2440 2440 344 1.12 181600 51240 12.5 1
4.064 3256 102 98 3 0 0</p>
<figure>
<img src="/OpenSeesRT/contrib/static/321.png"
title="Fig 2. Comparison between wood-sheathed CFS SWP experimental and numerical results specimens No 32"
width="400"
alt="Fig 2. Comparison between wood-sheathed CFS SWP experimental and numerical results specimens No 32" />
<figcaption aria-hidden="true">Fig 2. Comparison between wood-sheathed
CFS SWP experimental and numerical results specimens No 32</figcaption>
</figure>
<p>uniaxialMaterial CFSWSWP 16 2440 610 344 1.12 181600 51240 12.5 1
4.064 3256 152 40 3 0 0</p>
<figure>
<img src="/OpenSeesRT/contrib/static/16.png"
title="Fig 3. Comparison between wood-sheathed CFS SWP experimental and numerical results specimens No 16"
width="400"
alt="Fig 3. Comparison between wood-sheathed CFS SWP experimental and numerical results specimens No 16" />
<figcaption aria-hidden="true">Fig 3. Comparison between wood-sheathed
CFS SWP experimental and numerical results specimens No 16</figcaption>
</figure>
<hr />
<h2>References</h2>
<p><a
href="http://www.sciencedirect.com/science/article/pii/S0263823115301026">Smail
Kechidi and Nouredine Bourahla, Deteriorating hysteresis model for
cold-formed steel shear wall panel based on its physical and mechanical
characteristics, Journal of Thin-Walled Structures (2016), pp.421-430.
<a
href="DOI:10.1016/j.tws.2015.09.022">DOI:10.1016/j.tws.2015.09.022</a>.</a></p>
<p>Smail Kechidi, Hysteresis model development for cold-formed steel
shear wall panel based on physical and mechanical characteristics,
Master Thesis, University of Blida 1, Algeria, 2014.</p>
<p>Smail Kechidi, Nouredine Bourahla, Deteriorating hysteresis model for
cold-formed steel shear wall panel based on physical and mechanical
characteristics, OpenSees Days Portugal 2014- Workshop on Multi-Hazard
Analysis of Structures using OpenSees, Porto 3-4, Portugal, July
2014.</p>
<p>L.N. Lowes, A. Altoontash, Modelling reinforced-concrete beam-column
joints subjected to cyclic loading, Journal of Structural Engineering,
129(12):1686-1697, 2003.</p>
<p>A.E. Branston, Y.C. Chen, F.A. Boudreault and C.A. Rogers, Testing of
light-gauge steel frame wood structural panel shear walls, Canadian
Journal of Civil Engineering, 33(9):561-572, 2006.</p>
<p>J. Martinez and L. Xu, Strength and stiffness determination of shear
wall panels in cold-formed steel framing, Thin-Walled Structures,
44(10):1084-1095, 2006.</p>
<hr />
<p>Code Developed by: <span style="color:blue"> Smail Kechidi and
Nouredine Bourahla, University of Blida 1, Algeria </span></p>
<p>Images Developed by: <span style="color:blue"> Smail Kechidi,
University of Blida 1, Algeria </span></p>
<hr />
<p>Authors contact:</p>
<p><strong>Smail Kechidi</strong>, PhD student at University of Blida 1,
Algeria, s_kechidi@univ-blida.dz, skechidi@yahoo.com</p>
<p><strong>Nouredine Bourahla</strong>, Professor at University of Blida
1, Algeria, nbourahla@univ-blida.dz</p>
