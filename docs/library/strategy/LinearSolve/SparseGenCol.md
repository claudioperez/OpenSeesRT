# SparseGenColLinSolver

```cpp
#include <system_of_eqn/linearSOE/SparseGen/SparseGenColLinSolver.h>

class SparseGenColLinSolver: public LinearSOESolver
```


`SparseGenColLinSolver` is an abstract class. The `SparseGenColLinSolver`
class provides access for each subclass to the SparseGenColLinSOE object
through the pointer *theSOE*, which is a protected pointer.

### Constructor

\
### Destructor

\
// Public Methods\

\

The integer *classTag* is passed to the LinearSOESolver classes
constructor.

\
Does nothing, provided so the subclasses destructor will be called.

\
Sets the link to the SparseGenColLinSOE object *theSOE*. This is the
object on which the solver will perform the numerical computations.
