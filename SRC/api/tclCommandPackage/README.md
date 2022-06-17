# OpenSeesRT

## Features

- **Idempotent**
- **Direct access to components** Developing OpenSees components
  is made easier by providing a direct interface to core object
  types, allowing calls constructs like `material.getStress(0.002, commit=False)`
  - UniaxialMaterial
  - FiberSection

- **Structured problem representation**
  - Allows introspective add-ons ([Torsion]())
  - Better visualization tools, independent of core source code.

- Promotes stability of OpenSees core.
- Supports more python versions (3.6+)

Additional minor features:
- Verbosity control
- new `with` Tcl command and Python constructs

## Paradigm

- no more "interpreters"; An analysis 

## User Changes

When `OpenSeesRT` is loaded as a Tcl library, there are a few minor
changes from the classic `OpenSees` interpreter:

### Streams
- `puts` command prints to `stdout` by default, whereas classic OpenSees
  writes to `stderr`.

- new `redirect` command

- Dropped:

    // extern void *OPS_WFSection2d(G3_Runtime*);
    // extern void *OPS_RCCircularSection(G3_Runtime*);
    // extern void *OPS_RCSection2d(G3_Runtime*);
    // extern void *OPS_RCTBeamSection2d(G3_Runtime*);
    // extern void *OPS_RCTunnelSection(G3_Runtime*);
    // extern void *OPS_TubeSection(G3_Runtime*);

## Developer Changes

- No more `OPS_GetInt(void)`; use your host's API

- `ModelBuilder` namespacing functionality
  - Eliminates random code in important places like `Domain::Print`


## Codebase changes

Files that are superseded

- socket.c
- utility/PeerNGA.cpp

## Cleaning & TODO

Remove dependence on

- utility/SimulationInformation.\*
- utility/StringContainer.\*

Remove TimeSeriesIntegrators from C++; handle in pre-processing?


  CC="clang" CXX="clang++" cmake -DCMAKE_CXX_INCLUDE_WHAT_YOU_USE=include-what-you-use ..


