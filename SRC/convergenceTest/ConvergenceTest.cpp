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
// Purpose: This file contains the class definition for ConvergenceTest,
// which is an abstract class. Objects of concrete subclasses can be used
// to test the convergence of an algorithm.
//
// Written: fmk
// Date: 09/98
// Revised:
//
#include <ConvergenceTest.h>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <sstream>
 

ConvergenceTest::ConvergenceTest(int clasTag)
:MovableObject(clasTag)
{

}

ConvergenceTest::~ConvergenceTest()
{

}

std::string
ConvergenceTest::pad(double x)
{
  std::ostringstream oss;
  oss << std::setw(11) << x;
  return oss.str();
}


std::string
ConvergenceTest::pad(int i)
{
  const int n = 5; // i < 100 ? 3 : 5;
  std::ostringstream oss;
  oss << std::setw(n) << i;
  return oss.str();
}
