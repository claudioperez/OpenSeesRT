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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-04-04 16:55:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/matrixutil/MatrixUtil.cpp,v $
 
#include <math.h>
                                                                        
#include <stdlib.h>
#include <Vector.h>
#include "cbdi.h"

void invertMatrix(int n, const Matrix &a, Matrix &b)
{
  a.Invert(b);
}

void getCBDIinfluenceMatrix(int nIntegrPts, const Matrix &xi_pt, double L, Matrix &ls)
{
   // setup Vandermode and CBDI influence matrices
   double xi;
   Matrix G(nIntegrPts, nIntegrPts); 
   Matrix Ginv(nIntegrPts, nIntegrPts);
   Matrix l(nIntegrPts, nIntegrPts);
   Matrix I(nIntegrPts,nIntegrPts);      // an identity matrix for matrix inverse

   for (int i = 1; i <= nIntegrPts; i++)
      for (int j = 1; j <= nIntegrPts; j++)
      {
         int i0 = i - 1;
         int j0 = j - 1;
         xi = xi_pt(i0,0);
         G(i0,j0) =  pow(xi,j-1);
         l(i0,j0) = (pow(xi,j+1)-xi)/(j*(j+1));
      }
   
   I.Zero();
   for (int i=0; i<nIntegrPts; i++)
     I(i,i) = 1.0;

   //invertMatrix(nIntegrPts, G, Ginv);
   if (G.Solve(I,Ginv) < 0)
     opserr << "LargeDispBeamCol3d::getCBDIinfluenceMatrix() - could not invert G\n";
      
   // ls = l * Ginv * (L*L);
   ls.addMatrixProduct(0.0, l, Ginv, L*L);
}

void getCBDIinfluenceMatrix(int nIntegrPts, double *pts, double L, Matrix &ls)
{
   // setup Vandermode and CBDI influence matrices
   double xi;
   Matrix G(nIntegrPts, nIntegrPts); 
   Matrix Ginv(nIntegrPts, nIntegrPts);
   Matrix l(nIntegrPts, nIntegrPts);
   Matrix I(nIntegrPts,nIntegrPts);      // an identity matrix for matrix inverse

   for (int i = 0; i < nIntegrPts; i++) {
     xi = pts[i];
     for (int j = 1; j <= nIntegrPts; j++) {
       int j0 = j - 1;
       G(i,j0) =  pow(xi,j-1);
       l(i,j0) = (pow(xi,j+1)-xi)/(j*(j+1));
     }
   }
   
   I.Zero();
   for (int i=0; i<nIntegrPts; i++)
     I(i,i) = 1.0;

   //invertMatrix(nIntegrPts, G, Ginv);
   if (G.Solve(I,Ginv) < 0)
     opserr << "getCBDIinfluenceMatrix() - could not invert G\n";
      
   // ls = l * Ginv * (L*L);
   ls.addMatrixProduct(0.0, l, Ginv, L*L);
}

void getCBDIinfluenceMatrix(int nPts, double *pts, int nIntegrPts, double *integrPts, double L, Matrix &ls)
{
   // setup Vandermode and CBDI influence matrices
   double xi;
   Matrix G(nIntegrPts, nIntegrPts); 
   Matrix Ginv(nIntegrPts, nIntegrPts);
   Matrix l(nPts, nIntegrPts);
   Matrix I(nIntegrPts,nIntegrPts);      // an identity matrix for matrix inverse

   // Loop over columns
   for (int j = 1; j <= nIntegrPts; j++) {
     int j0 = j - 1;
     for (int i = 0; i < nIntegrPts; i++) {
       xi = integrPts[i];
       G(i,j0) =  pow(xi,j-1);
     }
     for (int i = 0; i < nPts; i++) {
       xi = pts[i];
       l(i,j0) = (pow(xi,j+1)-xi)/(j*(j+1));
     }
   }
   
   I.Zero();
   for (int i=0; i<nIntegrPts; i++)
     I(i,i) = 1.0;

   //invertMatrix(nIntegrPts, G, Ginv);
   if (G.Solve(I,Ginv) < 0)
     opserr << "getCBDIinfluenceMatrix() - could not invert G\n";
      
   // ls = l * Ginv * (L*L);
   ls.addMatrixProduct(0.0, l, Ginv, L*L);
}
