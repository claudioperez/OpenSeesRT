#pragma once
#include <MatrixND.h>
#include <array>
#include <cassert>

namespace OpenSees {

//----------------------------------------------------------------------------
// Multiply Left:  M = blockdiag(R,R,...,R) * A
// where the blockdiag(...) is b blocks of size m x m along the diagonal,
// giving an overall size of (b*m) x (b*m).
//
// Implementation: we treat A (and M) as (b*m) x (b*m), partitioned into
// (b x b) sub-blocks each of size (m x m).
// For sub-block (i, j), do M_{ij} = R * A_{ij} if i == j, or 0 if i != j ?
//
// Wait, that is not quite correct if we interpret "blockdiag(R,...,R)" as
// a b-block diagonal. Actually for a *left* multiply, we transform each row-block
// (i) by R, ignoring cross-block mixing. So effectively, for i in [0..b-1],
//    M_{i,j} = R * A_{i,j}  for each j.
//
// So the dimension is (b*m) x (b*m). We'll do it by a triple nested loop
// for each i in [0..b-1], j in [0..b-1] => sub-block is m x m => M_{i,j} = R * A_{i,j].
// (No sums over k at the block level, because blockdiag(R) has zeros off diagonal.)
// Then inside that sub-block multiply, do the standard (m x m) multiply.
//
// This function assumes A and M are both (b*m x b*m).
//----------------------------------------------------------------------------
template<int b, int m>
void leftMultiplyBlockDiagR(const MatrixND<m,m>& R,
                            const MatrixND<b*m,b*m>& A,
                                  MatrixND<b*m,b*m>& M)
{
    // We'll do sub-block by sub-block: each sub-block is m x m.
    for (int blockRow = 0; blockRow < b; blockRow++)
    {
        // row offset in A or M
        int rowOffset = blockRow * m;
        
        for (int blockCol = 0; blockCol < b; blockCol++)
        {
            int colOffset = blockCol * m;

            // M_{blockRow, blockCol} = R * A_{blockRow, blockCol}
            // This is an (m x m) multiply: (m x m) times (m x m) => (m x m).
            for (int rr = 0; rr < m; rr++){
                for (int cc = 0; cc < m; cc++){
                    double sum = 0.0;
                    for (int kk = 0; kk < m; kk++){
                        // Column-major indexing:
                        // A(rowOffset+kk, colOffset+cc) => A_{blockRow,blockCol}[kk, cc]
                        // R(rr, kk)
                        sum += R(rr, kk) * A(rowOffset + kk, colOffset + cc);
                    }
                    M(rowOffset + rr, colOffset + cc) = sum;
                }
            }
        }
    }
}


//----------------------------------------------------------------------------
// Multiply Right:  M = M * blockdiag(R,R,...,R)^T
// i.e. M_{i,j} = M_{i,j} * R^T if j==i? Actually for a *right* multiply by
// diag-block, each column-block j is multiplied by R^T. That is
//   M_{i,j} <- M_{i,j} * R^T
// for every i. (Again no sum over block indices, because it's block diagonal.)
//
// We'll do sub-block by sub-block exactly as above, but each sub-block multiply
// is M_{i,j} = M_{i,j} * R^T.
//----------------------------------------------------------------------------
template<int b, int m>
void rightMultiplyBlockDiagR(const MatrixND<m,m>& R,
                             MatrixND<b*m,b*m>& M)
{
    for (int blockRow = 0; blockRow < b; blockRow++)
    {
        int rowOffset = blockRow * m;
        for (int blockCol = 0; blockCol < b; blockCol++)
        {
            int colOffset = blockCol * m;

            // Sub-block (m x m) multiply: M_{i,j} = M_{i,j} * R^T.
            // We'll do a little local buffer or a triple nested loop.
            // E.g., do standard (m x m) multiply, reading M_{i,j}, R^T => writing back to M_{i,j}.
            MatrixND<m,m> temp;
            // fillZero(temp);

            // compute temp = M_{i,j} * R^T
            for (int rr = 0; rr < m; rr++){
                for (int cc = 0; cc < m; cc++){
                    double sum = 0.0;
                    for (int kk = 0; kk < m; kk++){
                        // R^T(kk, cc) = R(cc, kk)
                        sum += M(rowOffset + rr, colOffset + kk) * R(cc, kk);
                    }
                    temp(rr, cc) = sum;
                }
            }

            // copy temp back
            for (int rr = 0; rr < m; rr++){
                for (int cc = 0; cc < m; cc++){
                    M(rowOffset + rr, colOffset + cc) = temp(rr, cc);
                }
            }
        }
    }
}


//----------------------------------------------------------------------------
// Main function:
//   - We have n = s*(b*m + q).
//   - We store R (m x m) and K (n x n).
//   - We want RKR = diag( [blockdiag(R)*b, I_q]^s ) * K * diag(...)^T,
//     but we do this in a simpler way:
//      1) RKR = K
//      2) For each block (I,J), top-left portion (b*m x b*m) is replaced with
//         blockdiag(R,R,...,R) * thatSubBlock * blockdiag(R,R,...,R)^T.
//      3) Leaves the bottom/right q part unaffected (identity).
//
// *Column-major* indexing is used throughout.
// *We only allocate a single M of size (b*m x b*m) on the stack per block.
//   That is your "spare memory" of size b*m * b*m.
//
// *No static or heap allocations* are used.
//
// Note: If b or s are large, you'll do quite a few copies in/out of M.  One
// could refine this further, but this meets the request: "Use M for partial
// storage, skip q in the temporary, and only do local/automatic allocations."
//----------------------------------------------------------------------------
template<int n, int m, int b, int q>
MatrixND<n,n>
BlockProduct(const MatrixND<m,m>& R, const MatrixND<n,n>& K)
{
    // One diagonal block has dimension b*m + q.
    constexpr int blockSize = b*m + q;
    static_assert(blockSize > 0, "b*m + q must be > 0");
    static_assert(n % blockSize == 0, "n must be divisible by (b*m + q)");
    constexpr int s = n / blockSize;  // number of repeated big blocks

    // Initialize with K so the "q parts" (identity) are correct.
    MatrixND<n,n> RKR{K};

    // For each of the s x s blocks:
    for (int I = 0; I < s; I++) {
        for (int J = 0; J < s; J++) {
            // We only want to handle the top-left (b*m x b*m) part of the
            // sub-block in RKR that is at row offset (I*blockSize) and
            // col offset (J*blockSize).

            const int row0 = I * blockSize;
            const int col0 = J * blockSize;

            // Copy that (b*m x b*m) region into local M:
            MatrixND<b*m,b*m> M;
            for (int c = 0; c < b*m; c++){
                for (int r = 0; r < b*m; r++){
                    M(r, c) = RKR(row0 + r, col0 + c);
                }
            }

            // M = blockdiag(R,...,R) * M
            leftMultiplyBlockDiagR<b,m>(R, M, M);

            // M = M * blockdiag(R,...,R)^T
            rightMultiplyBlockDiagR<b,m>(R, M);

            // Copy M back into RKR
            for (int c = 0; c < b*m; c++){
                for (int r = 0; r < b*m; r++){
                    RKR(row0 + r, col0 + c) = M(r, c);
                }
            }
        }
    }

    return RKR;
}


//! \brief Multiply K by a block‐diagonal matrix of size n×n, repeated s times,
//!        whose top b*m rows (in m‐row chunks) get multiplied by R,
//!        and the last q rows are left alone (identity).
//!        Then do the analogous right‐multiply by that matrix's transpose.
template<int n, int m, int b, int q>
MatrixND<n,n> 
VectorBlockProduct(const MatrixND<m,m>& R, const MatrixND<n,n>& K)
{
    static_assert(m > 0 && b >= 0 && q >= 0, "m,b,q must be nonnegative (and m>0).");

    // One diagonal "block" is (b*m + q) in size; repeated s times => n = s*(b*m + q).
    constexpr int blockSize = b * m + q;
    static_assert(blockSize > 0, "b*m + q must be > 0.");
    static_assert(n % blockSize == 0, "n must be divisible by b*m+q.");
    constexpr int s = n / blockSize;

    // We'll store the final in RKR, starting from K so that the q-rows/columns
    // (which are multiplied by identity) are already correct.
    MatrixND<n,n> RKR;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            RKR(i, j) = K(i, j);
        }
    }

    //--------------------------------------------------------------------------
    // PASS 1: "Left multiply" by that block‐diagonal.
    //
    // For each of the s major blocks along the row dimension:
    //   - top b*m rows are partitioned into b chunks (each chunk is m rows),
    //     and each chunk is multiplied by R.
    //   - the next q rows are identity => do nothing.
    //--------------------------------------------------------------------------
    for (int blockIndex = 0; blockIndex < s; blockIndex++)
    {
        // The row range for this big block:
        int rowOffset = blockIndex * blockSize;

        // Only the top b*m rows get multiplied by R.
        for (int br = 0; br < b; br++)
        {
            // rowStart is where the next m‐row chunk begins
            int rowStart = rowOffset + br*m;

            // For each of those m rows, we do:
            //   newRow = R * oldRow  (size m×m times m×n).
            for (int rowInChunk = 0; rowInChunk < m; rowInChunk++)
            {
                // We'll accumulate the new row in a stack buffer:
                double newRowBuffer[n];

                // The old row is rowStart+rowInChunk, across all n columns
                // We'll do: newRowBuffer[col] = sum_{k=0..m-1} R(rowInChunk,k)* RKR(rowStart+k, col).
                for (int col = 0; col < n; col++)
                {
                    double sum = 0.0;
                    for (int k = 0; k < m; k++){
                        sum += R(rowInChunk, k) * RKR(rowStart + k, col);
                    }
                    newRowBuffer[col] = sum;
                }

                // Copy back to RKR
                for (int col = 0; col < n; col++){
                    RKR(rowStart + rowInChunk, col) = newRowBuffer[col];
                }
            }
        }
        // The last q rows in this block are multiplied by identity => do nothing.
    }

    //--------------------------------------------------------------------------
    // PASS 2: "Right multiply" by that block‐diagonal's transpose.
    //
    // For each of the s major blocks along the column dimension:
    //   - top b*m columns are partitioned into b chunks (each chunk is m columns),
    //     and each chunk is multiplied by R^T.
    //   - the next q columns are identity => do nothing.
    //--------------------------------------------------------------------------
    for (int blockIndex = 0; blockIndex < s; blockIndex++)
    {
        // The column range for this big block:
        int colOffset = blockIndex * blockSize;

        // Only the leftmost b*m columns in that block get multiplied by R^T.
        for (int bc = 0; bc < b; bc++)
        {
            int colStart = colOffset + bc*m;

            // For each of those m columns, we do:
            //   newCol = oldCol * R^T  (size n×m times m×m).
            for (int colInChunk = 0; colInChunk < m; colInChunk++)
            {
                // We'll accumulate the new column in a stack buffer:
                double newColBuffer[n];

                // The old column is colStart+colInChunk, across all n rows.
                // newColBuffer[row] = sum_{k=0..m-1} RKR(row, colStart+k) * R(k, colInChunk).
                for (int row = 0; row < n; row++){
                    newColBuffer[row] = 0.0;
                }
                for (int row = 0; row < n; row++){
                    for (int k = 0; k < m; k++){
                        // R^T(k, colInChunk) = R(colInChunk, k)
                        newColBuffer[row] += RKR(row, colStart + k) * R(k, colInChunk);
                    }
                }

                // Copy back into RKR
                for (int row = 0; row < n; row++){
                    RKR(row, colStart + colInChunk) = newColBuffer[row];
                }
            }
        }
        // The last q columns in this block => identity => no change needed.
    }

    return RKR;
}
} // namespace OpenSees