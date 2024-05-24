# ProfileSPDLinDirectThread

UNDER CONSTRUCTION



```cpp
#include <system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectThreadSolver.h>


class ProfileSPDLinDirectThreadSolver: public LinearSOESolver
```


A `ProfileSPDLinDirectThreadSolver` object can be constructed to solve a
ProfileSPDLinSOE object. It does this in parallel using threads by
direct means, using the $LDL^t$ variation of the Cholesky factorization.
The matrx $A$ is factored one row block at a time using a left-looking
approach. Within a row block the factorization is performed by $NP$
threads. No BLAS or LAPACK routines are called for the factorization or
subsequent substitution.

Constructor\

\
Destructor\

\
Public Methods\

\

\

\
A unique class tag (defined in  `<classTags.h>`) is passed to the
ProfileSPDLinSolver constructor.

\
Does nothing.

\
The solver first copies the B vector into X. FILL IN The solve process
changes $A$ and $X$.
