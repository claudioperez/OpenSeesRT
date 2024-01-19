# Cold-Formed Steel Steel Sheathed Shear Wall Panel examples

<p>In order to check the accuracy of the proposed models, results from
quasi-static non-linear analyses of Steel-Sheathed Cold-Formed Steel
Shear Wall Panel have been carried out using CFSSSWP uniaxialMaterial
model implemented into OpenSees version 2.4.6. For this purpose,
specimens 1C-b and 3C-a tested by Balh (2010) were selected from the
literature and analyzed under similar loading conditions. This set of
specimens covers a wide range of variation in physical and mechanical
characteristics such as: spacing, number, shear strength, diameter of
screw fasteners (sc, nc, Vs, and ds); wall aspect ratio (H/W); frame
thickness, chord stud moment of inertia, chord stud cross section area,
interior stud moment of inertia, yield and tensile strengths of steel
frame (tf, Af, Ife, Ifi, fyf, and fuf); sheathing thickness, yield and
tensile strengths of sheathing (type, ts, Fys, and Fus); as well as the
anchor bolt diameter of the HoldDown system (dt).</p>
<p>Each SWP was subjected to the Consortium of Universities for Research
in Earthquake Engineering (CUREE) loading protocol.</p>
<p><strong>TCL FILES:</strong></p>
<p><a href="Media:1Cbn.tcl" title="wikilink">1Cbn.tcl</a></p>
<p><a href="Media:3Can.tcl" title="wikilink">3Can.tcl</a></p>
<p><a href="Media:procRCycDAns.tcl"
title="wikilink">procRCycDAns.tcl</a></p>
<p>uniaxialMaterial CFSSSWP 1 2440 1220 346 396 1.14 436.22 395 300 0.46
1 4.166 1560 150 22.2 0 0</p>
<figure>
<img src="/OpenSeesRT/contrib/static/C1b.png"
title="Fig 1. Comparison between steel sheathed CFS SWP experimental and numerical results specimens No 1C-b"
width="400"
alt="Fig 1. Comparison between steel sheathed CFS SWP experimental and numerical results specimens No 1C-b" />
<figcaption aria-hidden="true">Fig 1. Comparison between steel sheathed
CFS SWP experimental and numerical results specimens No
1C-b</figcaption>
</figure>
<p>uniaxialMaterial CFSSSWP 2 2440 1220 391 342 0.87 436.22 395 300 0.46
1 4.166 1560 150 22.2 0 0</p>
<figure>
<img src="/OpenSeesRT/contrib/static/C3a.png"
title="Fig 2. Comparison between steel sheathed CFS SWP experimental and numerical results specimens No 3C-a"
width="400"
alt="Fig 2. Comparison between steel sheathed CFS SWP experimental and numerical results specimens No 3C-a" />
<figcaption aria-hidden="true">Fig 2. Comparison between steel sheathed
CFS SWP experimental and numerical results specimens No
3C-a</figcaption>
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
<p>Smail Kechidi and N Bourahla, Deteriorating hysteresis model for
cold-formed steel shear wall panel based on physical and mechanical
characteristics, OpenSees Days Portugal 2014- OPD 2014, 3-4 July 2014,
Porto, Portugal.</p>
<p>L.N. Lowes, A. Altoontash, Modelling reinforced-concrete beam-column
joints subjected to cyclic loading, Journal of Structural Engineering,
129(12):1686-1697, 2003.</p>
<p>Yanagi N, Yu C. Effective strip method for the design of cold-formed
steel framed shear wall with steel sheet sheathing. Journal of
Structural Engineering, ASCE 2014; 140(4).</p>
<p>Nisreen Balh, Development of seismic design provisions for steel
sheathed shear walls, Master Thesis, McGill University, Canada,
2010.</p>
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

