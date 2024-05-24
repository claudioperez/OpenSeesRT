# BandSPDLinLapackSolver

```tcl
system BandSPD
```


This command is used to construct a BandSPDSOE linear system of equation
object. As the name implies, this class is used for symmetric positive
definite matrix systems which have a banded profile. The matrix is
stored as shown below in a one dimensional array of size equal to the
(`bandwidth`/2) times the number of unknowns. When a solution is required,
the Lapack routines `DPBSV` and `DPBTRS` are used. 

## Theory

An *n*×*n* matrix *A*=(*a*~*i,j*~) is a **symmmetric banded matrix** if
all matrix elements are zero outside a diagonally bordered band whose
range is determined by constants *k*:

$$
a_{i,j}=0 \quad\texttt{if}\quad j<i-k \quad\texttt{ or }\quad j>i+k; \quad k \ge 0.$$

$$a_{i,j} = a_{j,i}$$

$$y^T A y  \ne 0 $$ 

for all non-zero vectors *y* with real entries
($y \in \mathbb{R}^n$),

The *bandwidth* of the matrix is *k* + *k* + 1.

For example, a symmetric 6-by-6 matrix with a right bandwidth of 2:

$$\begin{bmatrix}
 A_{11} & A_{12} & A_{13} &   0  & \cdots & 0 \\
      & A_{22} & A_{23} & A_{24} & \ddots & \vdots \\
      &        & A_{33} & A_{34} & A_{35} & 0 \\
      &        &        & A_{44} & A_{45} & A_{46} \\
      & sym    &        &        & A_{55} & A_{56} \\
      &        &        &        &        & A_{66}
\end{bmatrix}.$$ 

This matrix is stored as the 6-by-3 matrix:

$$\begin{bmatrix}
 A_{11} & A_{12} & A_{13} \\
 A_{22} & A_{23} & A_{24} \\
 A_{33} & A_{34} & A_{35} \\
 A_{44} & A_{45} & A_{46} \\
 A_{55} & A_{56} & 0 \\
 A_{66} & 0 & 0
\end{bmatrix}.
$$


## C++ Interface



A `BandSPDLinLapackSolver` object can be constructed to solve a
`BandSPDLinSOE` object. It obtains the solution by making calls on the the
LAPACK library. The class is defined to be a friend of the [`BandSPDLinSOE`](BandSPDLinSOE) class.




### Public Methods


A unique class tag (defined in  `<classTags.h>`) is passed to the
BandSPDLinSolver constructor.

\
Does nothing.

\
The solver first copies the B vector into X and then solves the
BandSPDLinSOE system by calling the LAPACK routines `dpbsv()`, if the
system is marked as not having been factored, and `dpbtrs()` if system
is marked as having been factored. If the solution is successfully
obtained, i.e. the LAPACK routines return $0$ in the INFO argument, it
marks the system has having been factored and returns $0$, otherwise it
prints a warning message and returns INFO. The solve process changes $A$
and $X$.

Does nothing but return $0$.

Does nothing but return $0$.

Does nothing but return $0$.

------------------------------------------------------------------------

Code Developed by: `<span style="color:blue">`{=html} fmk
`</span>`{=html}

