/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
#ifndef MatrixUtil_h
#define MatrixUtil_h

class Matrix;
void   invertMatrix(int n, const Matrix &a, Matrix &b);
void   getCBDIinfluenceMatrix(int nIntegrPts, const Matrix &xi_pt, double L, Matrix &ls);
void   getCBDIinfluenceMatrix(int nIntegrPts, double *pts, double L, Matrix &ls);
void   getCBDIinfluenceMatrix(int npts, double *pts, int nIntegrPts, double *ipts, double L, Matrix &ls);

#endif
