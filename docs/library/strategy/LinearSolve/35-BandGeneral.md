# BandGeneral SOE

<p>This command is used to construct a BandGeneralSOE linear system of
equation object. As the name implies, this class is used for matrix
systems which have a banded profile. The matrix is stored as shown below
in a 1dimensional array of size equal to the bandwidth times the number
of unknowns. When a solution is required, the Lapack routines DGBSV and
SGBTRS are used. The following command is used to construct such a
system:</p>

```tcl
system BandGeneral < -det >
```

<hr />

## Theory

An $n\times n$ matrix
$A=(a_{ij})$ is a <strong>band matrix</strong> if all matrix elements are zero outside a
diagonally bordered band whose range is determined by constants
<em>k</em><sub>1</sub> and
<em>k</em><sub>2</sub>:</p>
<dl>
<dt></dt>
<dd>

$$a_{i,j}=0 \quad\texttt{if}\quad j<i-k_1 \quad\texttt{ or
}\quad j>i+k_2; \quad k_1, k_2 \ge 0.$$

</dd>
</dl>
<p>The quantities <em>k</em><sub>1</sub> and
<em>k</em><sub>2</sub> are the <em>left</em> and
<em>right</em> <em>half-bandwidth</em>, respectively. The
<em>bandwidth</em> of the matrix is
<em>k</em><sub>1</sub>&nbsp;+&nbsp;<em>k</em><sub>2</sub>&nbsp;+&nbsp;1
(in other words, the smallest number of adjacent diagonals to which the
non-zero elements are confined).</p>
<p>and matrices are usually stored by storing the diagonals in the band;
the rest is implicitly zero.</p>

<p>For example, 6-by-6 a matrix with bandwidth 3:</p>

$$
\begin{bmatrix} B_{11} & B_{12} & 0 & \cdots & \cdots
& 0 \\ B_{21} & B_{22} & B_{23} & \ddots & \ddots
& \vdots \\ 0 & B_{32} & B_{33} & B_{34} & \ddots
& \vdots \\ \vdots & \ddots & B_{43} & B_{44} &
B_{45} & 0 \\ \vdots & \ddots & \ddots & B_{54} &
B_{55} & B_{56} \\ 0 & \cdots & \cdots & 0 & B_{65}
& B_{66} \end{bmatrix} 
$$

is stored as the 6-by-3 matrix

$$
\begin{bmatrix} 0 & B_{11} & B_{12}\\ B_{21} & B_{22}
& B_{23} \\ B_{32} & B_{33} & B_{34} \\ B_{43} & B_{44}
& B_{45} \\ B_{54} & B_{55} & B_{56} \\ B_{65} & B_{66}
& 0 \end{bmatrix}.
$$


# BandGenLinLapackSolver

This command is used to construct a BandGeneralSOE linear system of equation
object. As the name implies, this class is used for matrix systems which have a
banded profile. The matrix is stored as shown below in a one-dimensional array
of size equal to the bandwidth times the number of unknowns. When a solution is
required, the Lapack routines DGBSV and SGBTRS  are used.

## Theory

An $n\times n$ matrix *A*=(*a*~*i,j*~) is a **band matrix** if all matrix
elements are zero outside a diagonally bordered band whose range is
determined by constants *k*~1~ and *k*~2~:

$$a_{i,j}=0 \quad\texttt{if}\quad j<i-k_1 \quad\texttt{ or }\quad j>i+k_2; \quad k_1, k_2 \ge 0.$$

The quantities *k*~1~ and *k*~2~ are the *left* and *right*
*half-bandwidth*, respectively. The *bandwidth* of the matrix is
*k*~1~ + *k*~2~ + 1 (in other words, the smallest number of adjacent
diagonals to which the non-zero elements are confined).

and matrices are usually stored by storing the diagonals in the band;
the rest is implicitly zero.

For example, 6-by-6 a matrix with bandwidth 3:

$$
\begin{bmatrix}
 B_{11} & B_{12} & 0      & \cdots & \cdots & 0 \\
 B_{21} & B_{22} & B_{23} & \ddots & \ddots & \vdots \\
  0     & B_{32} & B_{33} & B_{34} & \ddots & \vdots \\
 \vdots & \ddots & B_{43} & B_{44} & B_{45} & 0 \\
 \vdots & \ddots & \ddots & B_{54} & B_{55} & B_{56} \\
 0      & \cdots & \cdots & 0      & B_{65} & B_{66}
\end{bmatrix}$$ 
is stored as the 6-by-3 matrix

$$\begin{bmatrix}
 0 & B_{11} & B_{12}\\
 B_{21} & B_{22} & B_{23} \\
 B_{32} & B_{33} & B_{34} \\
 B_{43} & B_{44} & B_{45} \\
 B_{54} & B_{55} & B_{56} \\
 B_{65} & B_{66} & 0
\end{bmatrix}.$$


## C++ Interface

```cpp
#include <system_of_eqn/linearSOE/bandGEN/BandGenLinLapackSolver.h>

class BandGenLinLapackSolver: public BandGenLinSolver
```

       MovableObject
       Solver
       LinearSOESolver


A `BandGenLinLapackSolver` object can be constructed to solve a
`BandGenLinSOE` object. It obtains the solution by making calls on the the
LAPACK library. The class is defined to be a friend of the BandGenLinSOE
class (see  `<BandGenLinSOE.h>`).


### Public Methods


A unique class tag (defined in  `<classTags.h>`) is passed to the
BandGenLinSolver constructor. Sets the size of *iPiv* to $0$, *iPiv*
being an integer array needed by the LAPACK routines.


Invokes delete on `iPiv` to free the memory allocated to store the
array.


The solver first copies the B vector into X and then solves the
BandGenLinSOE system by calling the LAPACK routines `dgbsv()`, if the
system is marked as not having been factored, and `dgbtrs()` if system
is marked as having been factored. If the solution is successfully
obtained, i.e. the LAPACK routines return $0$ in the INFO argument, it
marks the system has having been factored and returns $0$, otherwise it
prints a warning message and returns INFO. The solve process changes $A$
and $X$.

Is used to construct a 1d integer array, *iPiv* that is needed by the
LAPACK solvers. It checks to see if current size of *iPiv* is large
enough, if not it deletes the cold and creates a larger array. Returns
$0$ if successful, prints a warning message and returns a $-1$ if not
enough memory is available for this new array.


------------------------------------------------------------------------

<p>Code developed by: <span style="color:blue"> fmk</span></p>

