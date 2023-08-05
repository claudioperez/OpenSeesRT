#ifndef cmx_h
#define cmx_h

#ifdef __cplusplus
extern "C" {
#endif

int cmx_inv3   (double *a, double *ainv, int*ok_flag);
int cmx_inv4   (double *a, double *ainv, int*ok_flag);
int cmx_inv5   (double *a, double *ainv, int*ok_flag);
int cmx_inv6   (double *a, double *ainv, int*ok_flag);
int cmx_inv6_v2(double *a, double *ainv, int*ok_flag);
int cmx_inv6_v3(double *a, double *ainv, int*ok_flag);

#ifdef __cplusplus
}
#endif

#endif // cmx_h
