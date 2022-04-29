# OpenSeesRT

## Feature

- Idempotent
- Verbosity control
- Promotes stability of OpenSees core.

## User Changes

When `OpenSeesRT` is loaded as a Tcl library, there are a few minor
changes from the classic `OpenSees` interpreter:

### Streams
- `puts` command prints to `stdout` by default, whereas classic OpenSees
  writes to `stderr`.

- new `redirect` command

### Misc
- new `with` command





## Developer Changes

- No more `OPS_GetInt`; use your host's API

- `ModelBuilder` namespacing functionality
  - Eliminates random code in important places like `Domain::Print`

- ``

