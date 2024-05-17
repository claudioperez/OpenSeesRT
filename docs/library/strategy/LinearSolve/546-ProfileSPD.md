# ProfileSPD SOE

This command is used to construct a `ProfileSPDSOE` linear system of
equation. As the name implies, this class is used for symmetric
positive definite matrix systems. The matrix is stored as shown below in
a 1 dimensional array with only those values below the first non-zero
row in any column being stored. This is sometimes also referred to as a
skyline storage scheme. The following command is used to construct such
a system:

```tcl
system ProfileSPD 
```

<hr />
<p>THEORY:</p>
<p>An <em>n</em>&amp;times;<em>n</em> matrix
<em>A</em>=(<em>a</em>&lt;sub&gt;<em>i,j</em> &lt;/sub&gt;) is a
<strong>symmmetric postive definite matrix</strong> if:</p>
<dl>
<dt></dt>
<dd>

$$a_{i,j} = a_{j,i}\,$$

$y^T A y != 0$ for all non-zero vectors $y$ with real entries 
($y \in \mathbb{R}^n$.


In the skyline or profile storage scheme only the entries below the
first no-zero row entry in any column are stored if storing by rows: The
reason for this is that as no reordering of the rows is required in
gaussian eleimination because the matrix is SPD, no non-zero entries
will ocur in the elimination process outside the area stored.

For example, a symmetric 6-by-6 matrix with a structura as shown
below:

$$
\begin{bmatrix} 
A_{11} & A_{12} & 0 & 0 & 0 \\ 
  &      A_{22} & A_{23} & 0 & A_{25} \\ 
  &      &        A_{33} & 0 & 0 \\ 
  &      &               & A_{44} & A_{45} \\ 
  & sym & & & A_{55} \end{bmatrix}.$$

<p>The matrix is stored as 1-d array</p>

$$
\begin{bmatrix} A_{11} & A_{12} & A_{22} & A_{23} &
A_{33} & A_{44} & A_{25} & 0 & A_{45} & A_{55}
\end{bmatrix}.
$$

<p>with a further array containing indices of diagonal elements:</p>

$$
\begin{bmatrix} 1 & 3 & 5 & 6 & 10 \end{bmatrix}.
$$

<hr />
<p>Code Developed by: <span style="color:blue"> fmk
</span></p>
