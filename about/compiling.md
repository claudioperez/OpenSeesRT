

First install build dependencies:

- `cmake`
- `tcl`


| Dependency  | Package              |
|:------------|:---------------------|
| LAPACK      | `liblapack-dev`      |
| BLAS        | `libblas-dev`        |
| Tcl\*       | `tcl-dev`            |

To compile, 

- create a directory called `build/`
  ```shell
  mkdir build
  ```

- Configure the compilers:
  ```
  cd build
  cmakd ..
  ```

- Change back to the root directory (the one that contains `build/`) and run:
  ```shell
  make -C build OpenSees -j4
  ```
  
  
