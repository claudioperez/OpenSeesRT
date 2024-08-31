# About

The `opensees` package currently serves three objectives:

1. Executing existing OpenSees models
2. Easing the creation of OpenSees models, and
3. Interacting with OpenSees objects.

(1) is provided by the `tcl` and `openseespy` submodules, as
well as the command line interface which is accessed as follows:
```shell
python -m opensees
```

Experimental support for (2) is implemented in the `library` and `emit` submodules.

Experimental support for (3) is implemented in the `tcl.ModelRuntime` class,
as well as the `OpenSeesPyRT.so` C++ extension library.

