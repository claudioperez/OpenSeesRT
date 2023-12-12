//
// from https://github.com/code-saturne/code_saturne/blob/master/src/base/cs_math.h#L1025
//
void
cmx_invSY3_in_place(double  a[3][3])
{
  double a00 = a[1][1]*a[2][2] - a[2][1]*a[1][2];
  double a01 = a[2][1]*a[0][2] - a[0][1]*a[2][2];
  double a02 = a[0][1]*a[1][2] - a[1][1]*a[0][2];
  double a11 = a[0][0]*a[2][2] - a[2][0]*a[0][2];
  double a12 = a[1][0]*a[0][2] - a[0][0]*a[1][2];
  double a22 = a[0][0]*a[1][1] - a[1][0]*a[0][1];

  double det_inv = 1. / (a[0][0]*a00 + a[1][0]*a01 + a[2][0]*a02);

  a[0][0] = a00 * det_inv;
  a[0][1] = a01 * det_inv;
  a[0][2] = a02 * det_inv;
  a[1][0] = a01 * det_inv;
  a[1][1] = a11 * det_inv;
  a[1][2] = a12 * det_inv;
  a[2][0] = a02 * det_inv;
  a[2][1] = a12 * det_inv;
  a[2][2] = a22 * det_inv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the inverse of a symmetric matrix using Cramer's rule.
 *
 * \remark Symmetric matrix coefficients are stored as follows:
 *         (s11, s22, s33, s12, s23, s13)
 *
 * \param[in]     s      symmetric matrix
 * \param[out]    sout   sout = 1/s1
 */
/*----------------------------------------------------------------------------*/

void
cmx_invSY3(const double  s[6],
           double        sout[restrict 6])
{
  double detinv;

  sout[0] = s[1]*s[2] - s[4]*s[4];
  sout[1] = s[0]*s[2] - s[5]*s[5];
  sout[2] = s[0]*s[1] - s[3]*s[3];
  sout[3] = s[4]*s[5] - s[3]*s[2];
  sout[4] = s[3]*s[5] - s[0]*s[4];
  sout[5] = s[3]*s[4] - s[1]*s[5];

  detinv = 1. / (s[0]*sout[0] + s[3]*sout[3] + s[5]*sout[5]);

  sout[0] *= detinv;
  sout[1] *= detinv;
  sout[2] *= detinv;
  sout[3] *= detinv;
  sout[4] *= detinv;
  sout[5] *= detinv;
}

