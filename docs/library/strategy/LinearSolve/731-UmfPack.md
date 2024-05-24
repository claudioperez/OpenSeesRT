# UmfPack SOE

This command is used to construct a sparse system of equations which
uses the <a
href="http://www.cise.ufl.edu/research/sparse/umfpack/">UmfPack</a>
solver. The following command is used to construct such a system:

```tcl
system UmfPack  < -lvalueFact $LVALUE >
```

(LVALUE*the number of nonzero entries) is the amount of additional
memory set aside for fill in during the matrix solution, by default the
LVALUE factor is 10. You only need to experiment with this if you get
error messages back about LVALUE being too small.

<hr />

## References

<ul>
<li>A column pre-ordering strategy for the unsymmetric-pattern
multifrontal method, T. A. Davis, ACM Transactions on Mathematical
Software, vol 30, no. 2, June 2004, pp. 165-195.</li>
<li>Algorithm 832: UMFPACK, an unsymmetric-pattern multifrontal method,
T. A. Davis, ACM Transactions on Mathematical Software, vol 30, no. 2,
June 2004, pp. 196-199.</li>
<li>A combined unifrontal/multifrontal method for unsymmetric sparse
matrices, T. A. Davis and I. S. Duff, ACM Transactions on Mathematical
Software, vol. 25, no. 1, pp. 1-19, March 1999.</li>
<li>An unsymmetric-pattern multifrontal method for sparse LU
factorization, T. A. Davis and I. S. Duff, SIAM Journal on Matrix
Analysis and Applications, vol 18, no. 1, pp. 140-158, Jan. 1997.</li>
</ul>
<hr />
<p>Code Developed by: <span style="color:blue"> fmk
</span></p>
