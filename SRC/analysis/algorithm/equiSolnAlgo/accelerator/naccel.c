// naccel.f -- translated by f2c (version 20200916).
// You must link the resulting object file with libf2c:
//      on Microsoft Windows system, link with libf2c.lib;
//      on Linux or Unix systems, link with .../path/to/libf2c.a -lm
//      or, if you install libf2c.a in a standard place, with -lf2c -lm
//      -- in that order, at the end of the command line, as in
//              cc *.o -lf2c -lm
//      Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,
//
//              http://www.netlib.org/f2c/libf2c.zip

#include <tgmath.h>
// #include "f2c.h"

// Subroutine
int naccel_(int *n, int *itr, int *mvec, 
        double *tol, double *u, double *f)
{
    // System generated locals
    int u_dim1, u_offset, i__1, i__2;
    double d__1;

    // Builtin functions
    // double sqrt(double);

    // Local variables
    static double c__[10], h__[121];        // was [11][11]
    static int i__, j, k;
    static double t;
    static int tmp, head, nvec, link[11], last, next, jptr, kptr, km1ptr;

// ***********************************************************************

//   NACCEL -- Newton iteration accelerator.

//   Argument  I/O/M/T  Description
//   --------  -------  -----------
//     N          I     Vector size.

//     ITR        I     Newton iteration count.  Initialization is done
//                      for ITR=1 and ITR is ignored thereafter.

//     MVEC       I     Maximum number of vectors to use in the
//                      acceleration algorithm.  May change from call
//                      to call but should be greater than 0 and no
//                      more than the internal parameter MAX (=10).

//     TOL        I     Tolerance for dropping vectors.  We drop the
//                      pair (z_k,w_k) if the sine of the angle between
//                      w_k and span{w_1, ..., w_(k-1)} is less than TOL.

//     U          T     Work space of length at least N*(2*MVEC+2).
//                      It should not be modified between calls.

//     F          M     Vector of length N.  On entry, it is the value
//                      of the function f at the current iterate.  It
//                      is overwritten with the accelerated correction.
//                      The unaccelerated correction would simply be f
//                      itself.

// ***********************************************************************
//     ==================================================================
//      First call: save f and the (unaccelerated) correction.
//     ==================================================================
    // Parameter adjustments
    --f;
    u_dim1 = *n;
    u_offset = 1 + u_dim1 * 3;
    u -= u_offset;

    // Function Body
    if (*itr == 1) {
        head = 1;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            u[j + ((head << 1) + 1) * u_dim1] = f[j];
            u[j + ((head << 1) + 2) * u_dim1] = f[j];
// L10:
        }
        link[0] = 0;
        nvec = 1;
//       Free storage linked list.
        next = 2;
        for (k = 2; k <= 10; ++k) {
// L20:
            link[k - 1] = k + 1;
        }
        link[10] = 0;
        return 0;
    }
//     ==================================================================
//      Compute w_1.
//     ==================================================================
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
// L100:
        u[j + ((head << 1) + 2) * u_dim1] -= f[j];
    }
    t = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
// L101:
// Computing 2nd power
        d__1 = u[j + ((head << 1) + 2) * u_dim1];
        t += d__1 * d__1;
    }
    t = 1. / sqrt(t);
//     Normalize w_1 and apply same factor to z_1.
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
        u[j + ((head << 1) + 1) * u_dim1] = t * u[j + ((head << 1) + 1) * 
                u_dim1];
        u[j + ((head << 1) + 2) * u_dim1] = t * u[j + ((head << 1) + 2) * 
                u_dim1];
// L110:
    }
//     Update H.
    kptr = link[head - 1];
    i__1 = nvec;
    for (k = 2; k <= i__1; ++k) {
        h__[k * 11 - 11] = 0.;
        i__2 = *n;
        for (j = 1; j <= i__2; ++j) {
// L121:
            h__[k * 11 - 11] += u[j + ((head << 1) + 2) * u_dim1] * u[j + ((
                    kptr << 1) + 2) * u_dim1];
        }
        kptr = link[kptr - 1];
// L120:
    }
//     ==================================================================
//      Compute the Choleski factorization of H.
//     ==================================================================
    k = 2;
    h__[0] = 1.;
L200:
    if (k > fmin(nvec,*mvec)) {
        goto L250;
    }
    i__1 = k - 1;
    for (j = 1; j <= i__1; ++j) {
        h__[k + j * 11 - 12] = h__[j + k * 11 - 12];
        i__2 = j - 1;
        for (i__ = 1; i__ <= i__2; ++i__) {
// L211:
            h__[k + j * 11 - 12] -= h__[k + i__ * 11 - 12] * h__[j + i__ * 11 
                    - 12];
        }
        h__[k + j * 11 - 12] /= h__[j + j * 11 - 12];
// L210:
    }
    h__[k + k * 11 - 12] = 1.;
    i__1 = k - 1;
    for (j = 1; j <= i__1; ++j) {
// L220:
// Computing 2nd power
        d__1 = h__[k + j * 11 - 12];
        h__[k + k * 11 - 12] -= d__1 * d__1;
    }
// Computing 2nd power
    d__1 = *tol;
    if (h__[k + k * 11 - 12] < d__1 * d__1) {
//       -----------------------------------------------
//        w_k is nearly in span{w_1, ..., w_(k-1)}
//       -----------------------------------------------
//       Remove w_k from linked list.
        km1ptr = head;
        i__1 = k - 1;
        for (j = 2; j <= i__1; ++j) {
// L230:
            km1ptr = link[km1ptr - 1];
        }
        kptr = link[km1ptr - 1];
        link[km1ptr - 1] = link[kptr - 1];
        --nvec;
//       Update free storage list.
        link[kptr - 1] = next;
        next = kptr;
//       Update H.
        i__1 = nvec;
        for (j = k; j <= i__1; ++j) {
            i__2 = k - 1;
            for (i__ = 1; i__ <= i__2; ++i__) {
// L241:
                h__[i__ + j * 11 - 12] = h__[i__ + (j + 1) * 11 - 12];
            }
            i__2 = j - 1;
            for (i__ = k; i__ <= i__2; ++i__) {
// L242:
                h__[i__ + j * 11 - 12] = h__[i__ + 1 + (j + 1) * 11 - 12];
            }
// L240:
        }
        goto L200;
    } else {
        h__[k + k * 11 - 12] = sqrt(h__[k + k * 11 - 12]);
        ++k;
        goto L200;
    }
//     ------------------------------
//      Retain at most MVEC vectors.
//     ------------------------------
L250:
    if (nvec > *mvec) {
//       truncate the linked list.
        last = head;
        i__1 = *mvec;
        for (j = 2; j <= i__1; ++j) {
// L260:
            last = link[last - 1];
        }
        tmp = link[last - 1];
        link[last - 1] = 0;
        last = tmp;
//       Update free storage list.
        i__1 = nvec;
        for (j = *mvec + 2; j <= i__1; ++j) {
// L270:
            last = link[last - 1];
        }
        link[last - 1] = next;
        next = tmp;
        nvec = *mvec;
    }
//     ==================================================================
//      Compute the projection of f onto {w_1, ... , w_nvec}.
//     ==================================================================
    jptr = head;
    i__1 = nvec;
    for (j = 1; j <= i__1; ++j) {
        c__[j - 1] = 0.;
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
// L301:
            c__[j - 1] += f[i__] * u[i__ + ((jptr << 1) + 2) * u_dim1];
        }
        jptr = link[jptr - 1];
        i__2 = j - 1;
        for (i__ = 1; i__ <= i__2; ++i__) {
// L310:
            c__[j - 1] -= h__[j + i__ * 11 - 12] * c__[i__ - 1];
        }
        c__[j - 1] /= h__[j + j * 11 - 12];
// L300:
    }
    for (j = nvec; j >= 1; --j) {
        i__1 = nvec;
        for (i__ = j + 1; i__ <= i__1; ++i__) {
// L321:
            c__[j - 1] -= h__[i__ + j * 11 - 12] * c__[i__ - 1];
        }
        c__[j - 1] /= h__[j + j * 11 - 12];
// L320:
    }
//     ==================================================================
//      Compute the accelerated correction.
//     ==================================================================
//     Save f for the next call.
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
// L410:
        u[j + ((next << 1) + 2) * u_dim1] = f[j];
    }
    kptr = head;
    i__1 = nvec;
    for (k = 1; k <= i__1; ++k) {
        i__2 = *n;
        for (j = 1; j <= i__2; ++j) {
// L421:
            f[j] = f[j] - c__[k - 1] * u[j + ((kptr << 1) + 2) * u_dim1] + 
                    c__[k - 1] * u[j + ((kptr << 1) + 1) * u_dim1];
        }
        kptr = link[kptr - 1];
// L420:
    }
//     Save the correction for the next call.
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
// L430:
        u[j + ((next << 1) + 1) * u_dim1] = f[j];
    }
//     ==================================================================
//      Shift the vectors to the right.
//     ==================================================================
    tmp = next;
    next = link[tmp - 1];
    link[tmp - 1] = head;
    head = tmp;
//     Update H.
    for (j = nvec; j >= 1; --j) {
        i__1 = j - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
// L500:
            h__[i__ + 1 + (j + 1) * 11 - 12] = h__[i__ + j * 11 - 12];
        }
    }
    ++nvec;
    return 0;
} // naccel_

