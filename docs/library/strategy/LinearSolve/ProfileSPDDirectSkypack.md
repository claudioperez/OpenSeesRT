# ProfileSPDDirectSkypack

```cpp
#include <system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectSkypackSolver.h>

class ProfileSPDLinDirectSkypackSolver: public LinearSOESolver
```


A `ProfileSPDLinDirectSkypackSolver` object can be constructed to solve a
`ProfileSPDLinSOE` object. It does this by direct means using the routines
supplied in the SKYPACK library, a library which uses the BLAS levels
1,2 and 3 for the factorization and substitution.

The routines in SKYPACK require a number of work areas: *int block\[3$$
*
and *double invD\[size$$
*. In addition, to allow the use of the BLAS 2
and 3, work areas *double rw\[mRows x mCols$$
*, *double tw\[mRows x
mRows$$
* and *int index\[max(mCols,mRows)$$
* are created.


A unique class tag (defined in  `<classTags.h>`) is passed to the
ProfileSPDLinSolver constructor. Sets *mCols* and *mRows* equal to $0$
and does not try and allocate any memory for the work arrays.

A unique class tag (defined in  `<classTags.h>`) is passed to the
ProfileSPDLinSolver constructor. Sets *mCols* and *mRows* and allocates
space in memory for the work arrays *rw*, *tw* and *index*. If not
enough memory is available in memory, *mCols* and *mRows* is set equal
to $0$ and an error message is printed.

\
Invokes delete on any work areas that have been constructed.

\
The solver first copies the B vector into X and then solves the
BandSPDLinSOE system. If the matrix has not been factored, the matrix is
first factored using the SKYPACK routine *skysf2()*, if *mCols* and
*mRows* equal $0$, or *skypf2()*. *skysf2()* is a routine which uses the
BLAS level 1 routines, *skypf2()* is a routine which uses BLAS levels 2
and 3. If *skypf2()* has been called, *invD* is set up. Once the matrix
has been factored, `skyss()` is called. If the solution is successfully
obtained, i.e. the `skyss()` routine returns $0$ in the INFO argument,
$0$ is returned, otherwise it prints a warning message and returns INFO.
The solve process changes $A$ and $X$.

Is responsible for setting the *block* information required by the
SKYPACK routines (block\[0$$
=1; block\[1$$
=size, block\[2$$
=1) and for
creating space for the *invD* work array. Returns $0$ if successful,
otherwise a warning message is printed and a $-1$ is returned.

Does nothing but return $0$.

Does nothing but return $0$.
