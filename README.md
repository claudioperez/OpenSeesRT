# `opensees`



<!--
<img align="center" src="https://raw.githubusercontent.com/BRACE2/OpenSeesRT/master/docs/figures/banner.png" width="150px" alt="OpenSees">
-->

<img align="left" src="https://raw.githubusercontent.com/claudioperez/sdof/master/docs/assets/peer-black-300.png" width="150px" alt="PEER Logo">


The fastest package for nonlinear finite element analysis in the west.

<br>


<br>

<div style="align:center">

<!--
[![Latest conda-forge version](https://img.shields.io/conda/vn/conda-forge/opensees?logo=conda-forge&style=for-the-badge)](https://anaconda.org/conda-forge/opensees)
-->

<!-- [![PyPI Downloads][pypi-v-image]][pypi-v-link] -->

[![Latest PyPI version](https://img.shields.io/pypi/v/opensees?logo=pypi&style=for-the-badge)](https://pypi.python.org/pypi/opensees)
[![](https://img.shields.io/conda/v/opensees/opensees?color=%23660505&style=for-the-badge)](https://anaconda.org/opensees/opensees)
[![PyPI Downloads](https://img.shields.io/pypi/dm/opensees?style=for-the-badge)](https://pypi.org/project/opensees)

</div>

<!-- 
-------------------------------------------------------------------- 
-->


# OpenSees

`opensees` is a Python package that provides an intuitive API for nonlinear
finite element analysis, implemented in C++. The library features
state-of-the-art finite element formulations and solution algorithms, including
mixed formulations for beams and solids, over 200 material models, and an
extensive collection of continuation algorithms to solve highly nonlinear
problems. 


The package may be used as a drop-in replacement for both OpenSees.exe and
OpenSeesPy (see *Getting Started* below), and generally provides a substantial performance boost.

## Features

- **Interactive Tasks**: Easily return stiffness, mass, and damping matrices as NumPy arrays and join meshes without duplicate nodes and constraints.
- **Intuitive API**: Easy-to-use Python interface for complex finite element analysis.
- **Advanced Formulations**: Mixed formulations for beams and solids.
- **Extensive Material Models**: Over 200 material models to choose from.
- **Continuation Algorithms**: Robust algorithms for solving highly nonlinear problems.

- **Performance** The `opensees` package uses the experimental 
  [`OpenSeesRT`](https://github.com/claudioperez/OpenSeesRT) 
  analysis kernel which
  eliminates reliance on global variables for state and memory management. 
  New template classes in [`matrix/`](https://github.com/STAIRLab/OpenSeesRT/tree/master/SRC/matrix/)
  allow for stack-allocated
  matrices and vectors and eliminate the need for static local variables 
  which is leveraged to furnish substantial performance improvements. 
  Switching Python scripts
  to use `opensees` typically results in a 4x to 5x performance boost.

- **Reliability** The core OpenSees runtime has been redesigned so that all program 
  state is encapsulated in user-instantiated classes,
  and global variables/singletons are avoided. 
  This eliminates several preexisting vulnerabilities to inadvertent state corruption.


<!-- 
- **Semantics** Unlike interfaces which rely on global state, this package can be used 
  with true library semantics. 
-->

Additional features include:

- Convert OpenSeesPy scripts into equivalent Tcl files that can be used
  for faster processing or serialization. Unlike most conversion utilities,
  this conversion is done *exactly* and does not rely on hand-rolled parsing.

- The package can be installed with `pip` for Python versions 3.7 - 3.12 on Linux, MacOS and
  Windows, but eigenvalue analysis is currently broken on Windows.

> [!NOTE]
> This package is independent of the [`openseespy`](https://pypi.org/project/openseespy)
> library, which is documented in the OpenSees [documentation](https://opensees.github.io/OpenSeesDocumentation)
> website. OpenSeesPy can be installed by running the following command:
>
> ```shell
> pip install opensees[py]
> ```


<p style="text-align: center;">
<b>This package is <i>experimental</i> and not yet intended for public use.</b>
</p>


### Getting Started

- To start a Tcl interpreter run the shell command:

  ```bash
  python -m opensees
  ```

- The `opensees` package exposes a compatibility layer that exactly reproduces
  the original *OpenSeesPy* functions, but does so without mandating a single
  global program state. To run OpenSeesPy scripts, just change the import:
  ```python
  import openseespy.opensees
  ```
  to
  ```python
  import opensees.openseesrt
  ```
  For true stateless modeling, the `Model` class should be used instead of the legacy
  `model` function; documentation is under development.


## Compiling

To compile the project see [help/compiling](https://github.com/claudioperez/opensees/blob/master/help/compiling.md)

<!-- Badge links -->

[pypi-d-image]: https://img.shields.io/pypi/dm/opensees.svg
[license-badge]: https://img.shields.io/pypi/l/opensees.svg
[pypi-d-link]: https://pypi.org/project/opensees
[pypi-v-image]: https://img.shields.io/pypi/v/opensees.svg
[pypi-v-link]: https://pypi.org/project/opensees


See also:

- [`osmg`](https://pypi.org/project/osmg) OpenSees Model Generator
- [`sees`](https://pypi.org/project/sees) Modern rendering library
- [`mdof`](https://pypi.org/project/mdof) Optimized system identification library
- [`sdof`](https://pypi.org/project/sdof) Optimized integration for single degree of freedom systems

For more projects by the STAIRlab, visit https://github.com/STAIRlab .

## Support

<table align="center" style="border: 0;">
 <tr style="background-color:rgba(0, 0, 0, 0);">
  <td style="background-color:rgba(0, 0, 0, 0);" colspan="3">
    <a>
    <img src="https://raw.githubusercontent.com/claudioperez/opensees/master/docs/figures/opensees.svg" 
         width="600" alt="OpenSeesRT Logo">
    </a>
  </td>
 </tr>

<tr>
  <td>
    <a href="https://peer.berkeley.edu">
    <img src="https://raw.githubusercontent.com/claudioperez/sdof/master/docs/assets/peer-black-300.png"
         alt="PEER Logo" width="200"/>
    </a>
  </td>

  <td>
    <a href="https://dot.ca.gov/">
    <img src="https://raw.githubusercontent.com/claudioperez/sdof/master/docs/assets/Caltrans.svg.png"
         alt="Caltrans Logo" width="200"/>
    </a>
  </td>

  <td>
    <a href="https://brace2.herokuapp.com">
    <img src="https://raw.githubusercontent.com/claudioperez/sdof/master/docs/assets/stairlab.svg"
         alt="STAIRlab Logo" width="200"/>
    </a>
  </td>
 
 </tr>
</table>

