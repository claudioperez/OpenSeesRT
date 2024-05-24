# FullGeneral

This command is used to construct a Full General linear system of
equation object. As the name implies, the class utilizes NO space saving
techniques to cut down on the amount of memory used. If the matrix is of
size, nxn, then storage for an nxn array is sought from memory when the
program runs. When a solution is required, the Lapack routines `DGESV` and
`DGETRS` are used. The following command is used to construct such a
system:

```tcl
system FullGeneral < -det >
```

<hr />

<p>NOTES:</p>
<p>This type of system should almost never be used! This is because it
requires a lot more memory than every other solver and takes more time
in the actal solving operation than any other solver. It is required if
the user is interested in looking at the global system matrix.</p>

# FullGenLinLapackSolver

```cpp
#include<~/system_of_eqn/linearSOE/fullGEN/FullGenLinLapackSolver.h>

class FullGenLinLapackSolver: public FullGenLinSolver
```

    MovableObject
    Solver
    LinearSOESolver


A `FullGenLinLapackSolver` object can be constructed to solve a
FullGenLinSOE object. It obtains the solution by making calls on the the
LAPACK library. The class is defined to be a friend of the FullGenLinSOE
class (see  [`FullGenLinSOE`](FullGenLinSOE)).

### Constructor


### Destructor


### Public Methods


A unique class tag (defined in  `<classTags.h>`) is passed to the
FullGenLinSolver constructor. Sets the size of *iPiv* to $0$, *iPiv*
being an integer array needed by the LAPACK routines.

\
Invokes delete on *iPiv* to free the memory it was allocated.

\
First copies $B$ into $X$ and then solves the FullGenLinSOE system it is
associated with (pointer kept by parent class) by calling the LAPACK
routines `dgesv()`, if the system is marked as not having been factored,
or `dgetrs()`, if system is marked as having been factored. If the
solution is successfully obtained, i.e. the LAPACK routines return $0$
in the INFO argument, it marks the system has having been factored and
returns $0$, otherwise it prints a warning message and returns INFO. The
solve process changes $A$ and $X$.

Is used to construct a 1d integer array, *iPiv* that is needed by the
LAPACK solvers. It checks to see if current size of `iPiv` is large
enough, if not it deletes the cold and creates a larger array. Returns
$0$ if successful, prints a warning message and returns a $-1$ if not
enough memory is available for this new array.


<p>Code Developed by: <span style="color:blue">fmk</span></p>
