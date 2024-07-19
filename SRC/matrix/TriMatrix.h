//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//

#ifndef _OPS_TriDiagonalMatrixF
#define _OPS_TriDiagonalMatrixF
class TriDiagonalMatrixF
{    
  public: 
    // Construct an NxN matrix.
    TriDiagonalMatrixF(int n);
    ~TriDiagonalMatrixF();
    int N();
    void SetMat(int row, int col, double value);
    double GetMat(int row, int col);
    double* Solve(double* d, int dLength);

  private:
    int length;
    // The values for the sub-diagonal. A[0] is never used.
    double *A;

    // The values for the main diagonal.
    double *B;

    // The values for the super-diagonal. C[C.Length-1] is never used.
    double *C;
};
#endif

