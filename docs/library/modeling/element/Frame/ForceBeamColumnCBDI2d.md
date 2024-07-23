# ForceBeamColumnCBDI2d

:::
```tcl
element ForceBeamColumnCBDI eleTag iNode jNode transfTag integrationTag
        <-iter $maxIter> <-mass $mass>
```
:::


  // $Revision: 1.13 $
  // $Date: 2009-02-05 16:28:20 $
  // $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ForceBeamColumnCBDI2d.h,v $

## References
#### State Determination Algorithm
- Neuenhofer, A. and F. C. Filippou (1997). "Evaluation of Nonlinear Frame Finite
  Element Models." Journal of Structural Engineering, 123(7):958-966.

- Spacone, E., V. Ciampi, and F. C. Filippou (1996). "Mixed Formulation of
  Nonlinear Beam Finite Element." Computers and Structures, 58(1):71-83.


#### Plastic Hinge Integration
- Scott, M. H. and G. L. Fenves (2006). "Plastic Hinge Integration Methods for
  Force-Based Beam-Column Elements." Journal of Structural Engineering,
  132(2):244-252.


#### Analytical Response Sensitivity (DDM)
- Scott, M. H., P. Franchin, G. L. Fenves, and F. C. Filippou (2004).
  "Response Sensitivity for Nonlinear Beam-Column Elements."
  Journal of Structural Engineering, 130(9):1281-1288.


#### Software Design
- Scott, M. H., G. L. Fenves, F. T. McKenna, and F. C. Filippou (2007).
  "Software Patterns for Nonlinear Beam-Column Models."
  Journal of Structural Engineering, Approved for publication, February 2007.


