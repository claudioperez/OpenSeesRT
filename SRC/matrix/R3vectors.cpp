/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Ed "C++" Love
//
// #include <R3vectors.h>
#include <Vector.h>
#include <Matrix.h> 
#include <math.h>

#define sign(a) ( (a)>0 ? 1:-1 )

Vector  LovelyCrossProduct( const Vector &v, const Vector &w )
{

  Vector cross(3) ;

  cross(0) = v(1)*w(2) - v(2)*w(1) ;
  
  cross(1) = v(2)*w(0) - v(0)*w(2) ;

  cross(2) = v(0)*w(1) - v(1)*w(0) ;
 
  return cross ;

}

