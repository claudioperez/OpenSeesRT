# SparseSYM

This command is used to construct a sparse symmetric system of
equations which uses a row-oriented solution method in the solution
phase. The following command is used to construct such a system:

```tcl
system SparseSYM
```

<hr />
<p>NOTES:</p>
<ul>
<li>versions upto and including 2.2.0 use SparseSPD instead of SparseSYM
as the option to the system command. For backward compatability this old
option continues to work.</li>
</ul>
<hr />

## References
<p>Kincho H. Law and David R. McKay, “A Parallel Row-Oriented Sparse
Solution Method for Finite Element Structural Analysis,” International
Journal for Numerical Methods in Engineering, 36:2895-2919, 1993.</p>
<hr />
<p>Code Developed by: <span style="color:blue"> Jun Peng, Stanford
and now GOOGLE </span></p>
