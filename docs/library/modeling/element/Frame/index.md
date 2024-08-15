
# Frame

```
element Frame $tag $nodes $section $geometry
```

| Element                   | Shear | $P-\delta$ | Warp |  |
|---------------------------|-------|------------|------|--|
| `PrismFrame`              | y/n   |   0-5      | 0/1/2 (Inter)
| [`ForceFrame`](./force/)  | y/n   |   0/1
| `ShearFrame`              | y/n   |   0/1
| `EulerFrame`              | y/n   |   0/1
| `ExactFrame`              |   y   | Exact
| `MixedFrame`              | y/n   |   0/1

<!--
| `HingeFrame`      |       |
| `Elastica1D`      |
-->


# Transform Library

```
Logarithm   : Init | Incr | None
Geometry    : Linear      | P-Delta  | Corotational | Geodesic
Torsion[7?] : Wagner | None
Logarithm   : Init | Incr | None
```

