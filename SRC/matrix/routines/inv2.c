
#include <math.h>
#include <stdbool.h>

int cmx_inv2(double *a, double *ainv, int *ok_flag__)
{

/* **************************************************************************************** */
/*  m22inv  -  compute the inverse of a 2x2 matrix. */

/*  a       = input 2x2 matrix to be inverted */
/*  ainv    = output 2x2 inverse of matrix a */
/*  ok_flag = (output) .true. if the input matrix could be inverted, and .false. if the input matrix is singular. */
/* **************************************************************************************** */

    static double cofactor[4]	/* was [2][2] */;
    static int i__, j;
    static double det, eps;


    /* Parameter adjustments */
    ainv -= 3;
    a -= 3;

    /* Function Body */
    eps = 1e-10;
    det = a[3] * a[6] - a[5] * a[4];
    if (fabs(det) <= eps) {
	*ok_flag__ = false;
	return 0;
    }
    cofactor[0] = a[6];
    cofactor[2] = -a[4];
    cofactor[1] = -a[5];
    cofactor[3] = a[3];
    for (i__ = 1; i__ <= 2; ++i__) {
	for (j = 1; j <= 2; ++j) {
	    ainv[j + (i__ << 1)] = cofactor[i__ + (j << 1) - 3] / det;
	}
    }
    *ok_flag__ = true;
    return 0;
} /* m22inv_ */

