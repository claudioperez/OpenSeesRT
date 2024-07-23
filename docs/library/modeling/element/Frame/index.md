
# Frame

| Element           | Shear | $P-\delta$ | Warp |  |
|-------------------|-------|------------|------|--|
| `PrismFrame`      | y/n   |   0-5      | 0/1/2 (Inter)
| `DisplFrame`      | y/n   |   0/1
| `ForceFrame`      | y/n   |   0/1
| `ExactFrame`      |   y   | Exact
| `MixedFrame`      | y/n   |   0/1
| `MixedFrame3d_IW` |  y    |            | Intra (I/IJ, S/L)
| `MixedFrame3d_EW` | yn    |            | Inter
| `ExactFrame3d_EW` |  y    |            | Inter 
| `HingeFrame`      |       |
| `Elastica1D`      |



# Transform Library

```
Logarithm   : Init | Incr | None
Geometry    : Linear      | P-Delta  | Corotational | Geodesic
Torsion[7?] : Wagner | None
Logarithm   : Init | Incr | None
```

# Section Library

```
FiberSection
ShearOpenSection_wEW
ShearCellSection_wEW
ShearRectSection_wEW
ShearOpenSection_wIW
ShearCellSection_wIW
ShearRectSection_wIW
```


```tcl
section Fiber -Cw -GA -GJ {
    patch
    layer
    fiber
    ...
}

section Shear -warp [intra|inter] {
    layer 
    ...
}
```
