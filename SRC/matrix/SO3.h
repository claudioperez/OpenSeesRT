#ifndef SO3_H
#define SO3_H

namespace SO3 {
    /**
    * Computes the Spin of the input vector V, and saves the result into the output matrix S.
    * Note: no check is made on the size of the input-output arguments.
    * @param V the input vector (assumed size: >= 3)
    * @param S the output matrix (assumed size: >= 3x3)
    */
    template< class TVec, class TMat>
    inline static void Spin(const TVec& V, TMat& S)
    {
        S(0, 0) = 0.00;        S(0, 1) = -V[2];       S(0, 2) = V[1];
        S(1, 0) = V[2];        S(1, 1) =  0.00;       S(1, 2) = -V[0];
        S(2, 0) = -V[1];       S(2, 1) =  V[0];       S(2, 2) = 0.00;
    }

    /**
    * Computes the Spin of the input vector V, and saves the result into the output matrix S,
    * at the specified row index.
    * Note: no check is made on the size of the input-output arguments.
    * @param V the input vector (assumed size: >= 3)
    * @param S the output matrix (assumed size: >= 3x3)
    * @param row_index the index of the first row in the output matrix where the spin has to be saved
    */
    template< class TVec, class TMat>
    inline static void Spin_AtRow(const TVec& V, TMat& S, size_t row_index)
    {
        size_t i0 = row_index;
        size_t i1 = 1 + row_index;
        size_t i2 = 2 + row_index;
        double v0 = V(i0);
        double v1 = V(i1);
        double v2 = V(i2);
        S(i0, 0) = 0.00;    S(i0, 1) = -v2;        S(i0, 2) = v1;
        S(i1, 0) = v2;        S(i1, 1) = 0.00;    S(i1, 2) = -v0;
        S(i2, 0) = -v1;        S(i2, 1) = v0;        S(i2, 2) = 0.00;
    }

    /**
    * Computes the Spin of the input vector V, from the specified index, and saves the result into the output matrix S,
    * at the specified row index.
    * Note: no check is made on the size of the input-output arguments.
    * @param V the input vector (assumed size: >= 3)
    * @param S the output matrix (assumed size: >= 3x3)
    * @param vector_index the index of the first component of the input vector to be used to compute the spin
    * @param row_index the index of the first row in the output matrix where the spin has to be saved
    */
    template< class TVec, class TMat>
    inline static void Spin_AtRow(const TVec& V, TMat& S, size_t vector_index, size_t matrix_row_index)
    {
        size_t i0 = matrix_row_index;
        size_t i1 = 1 + matrix_row_index;
        size_t i2 = 2 + matrix_row_index;
        double v0 = V(vector_index);
        double v1 = V(vector_index + 1);
        double v2 = V(vector_index + 2);
        S(i0, 0) = 0.00;    S(i0, 1) = -v2;        S(i0, 2) = v1;
        S(i1, 0) = v2;        S(i1, 1) = 0.00;    S(i1, 2) = -v0;
        S(i2, 0) = -v1;        S(i2, 1) = v0;        S(i2, 2) = 0.00;
    }

    /**
    * Computes the Spin of the input vector V, and saves the result into the output matrix S.
    * This version uses a multiplier for the output values.
    * Note: no check is made on the size of the input-output arguments.
    * @param V the input vector (assumed size: >= 3)
    * @param S the output matrix (assumed size: >= 3x3)
    * @param mult the multiplier for the output values
    */
    template< class TVec, class TMat>
    inline static void Spin(const TVec& V, TMat& S, double mult)
    {
        S(0, 0) = 0.00;            S(0, 1) = -mult * V(2);    S(0, 2) = mult * V(1);
        S(1, 0) = mult * V(2);    S(1, 1) = 0.00;            S(1, 2) = -mult * V(0);
        S(2, 0) = -mult * V(1);    S(2, 1) = mult * V(0);    S(2, 2) = 0.00;
    }

    /**
    * Computes the Spin of the input vector V, and saves the result into the output matrix S,
    * at the specified row index.
    * This version uses a multiplier for the output values.
    * Note: no check is made on the size of the input-output arguments.
    * @param V the input vector (assumed size: >= 3)
    * @param S the output matrix (assumed size: >= 3x3)
    * @param mult the multiplier for the output values
    * @param row_index the index of the first row in the output matrix where the spin has to be saved
    */
    template< class TVec, class TMat>
    inline static void Spin_AtRow(const TVec& V, TMat& S, double mult, size_t row_index)
    {
        size_t i0 = row_index;
        size_t i1 = 1 + row_index;
        size_t i2 = 2 + row_index;
        double v0 = mult * V(i0);
        double v1 = mult * V(i1);
        double v2 = mult * V(i2);
        S(i0, 0) = 0.00;    S(i0, 1) = -v2;        S(i0, 2) = v1;
        S(i1, 0) = v2;        S(i1, 1) = 0.00;    S(i1, 2) = -v0;
        S(i2, 0) = -v1;        S(i2, 1) = v0;        S(i2, 2) = 0.00;
    }

    /**
    * Computes the Spin of the input vector V, from the specified index, and saves the result into the output matrix S,
    * at the specified row index.
    * This version uses a multiplier for the output values.
    * Note: no check is made on the size of the input-output arguments.
    * @param V the input vector (assumed size: >= 3)
    * @param S the output matrix (assumed size: >= 3x3)
    * @param mult the multiplier for the output values
    * @param vector_index the index of the first component of the input vector to be used to compute the spin
    * @param row_index the index of the first row in the output matrix where the spin has to be saved
    */
    template< class TVec, class TMat>
    inline static void Spin_AtRow(const TVec& V, TMat& S, double mult, size_t vector_index, size_t matrix_row_index)
    {
        size_t i0 = matrix_row_index;
        size_t i1 = 1 + matrix_row_index;
        size_t i2 = 2 + matrix_row_index;
        double v0 = mult * V(vector_index);
        double v1 = mult * V(vector_index + 1);
        double v2 = mult * V(vector_index + 2);
        S(i0, 0) = 0.00;    S(i0, 1) = -v2;        S(i0, 2) = v1;
        S(i1, 0) = v2;      S(i1, 1) = 0.0;        S(i1, 2) = -v0;
        S(i2, 0) = -v1;     S(i2, 1) =  v0;        S(i2, 2) = 0.00;
    }

};
#endif // SO3_H
