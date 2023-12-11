/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Ed "C++" Love
//
#include <R3vectors.h>
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


