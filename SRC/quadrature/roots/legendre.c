/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Gauss-Legendre integration function, gauleg, from "Numerical Recipes in C"
// (Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and
// B.P. Flannery
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// EPS is the relative precision.
#define EPS 3.0e-11

void legendre(double x1, double x2, double x[], double w[], int n);

int main(int argc, char** argv)
{
  double pt[100],
         wt[100];

  int n = atoi(argv[1]);

  legendre(-1, 1, pt, wt, n);

  for (int i=0; i < n; i++)
    printf("%lf\t%lf\n", pt[i], wt[i]);
}

void legendre(double x1, double x2, double x[], double w[], int n)
/*******************************************************************************
 Given the lower and upper limits of integration x1 and x2, and given n, this
 routine returns arrays x[n] and w[n] of length n, containing the abscissas
 and weights of the Gauss-Legendre n-point quadrature formula.
*******************************************************************************/
{
  // The roots are symmetric, so we only find half of them.
  int    m  = (n+1)/2; 
  double xm = 0.5*(x2+x1);
  double xl = 0.5*(x2-x1);

  // Loop over the desired roots. Note iteration begins at i=1
  for (int i=1; i<=m; i++) { 
    // Initialize z with approximation to the ith root
    double z1, pp;
    double z = cos(3.141592654*(i-0.25)/(n+0.5));

    // Commence Newton iterations
    do {
      // Compute Legendre polynomial recurrence relation.
      double p1 = 1.0;
      double p2 = 0.0;
      for (int j=1; j<=n; j++) { 
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j;
      }
      // p1 is now the desired Legendre polynomial. We next compute
      // its derivative in `pp` by a standard relation involving
      // p2, the polynomial of one lower order.
      pp = n*(z*p1 - p2)/(z*z - 1.0);
      z1 = z;
      z  = z1 - p1/pp; // Newton update.
    } while (fabs(z-z1) > EPS);

    x[i-1] = xm - xl*z;      // Scale the root to the desired interval, 
    x[n-i] = xm + xl*z;      // and put in its symmetric counterpart.
    w[i-1] = 2.0*xl/((1.0 - z*z)*pp*pp); // Compute the weight
    w[n-i] = w[i-1];                     // and its symmetric counterpart.
  }
}
