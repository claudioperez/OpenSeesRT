#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "cmx.h"

#define DIM 5
#define sizeDoubleWork 36

# define  DGESV  dgesv_   
# define  DGETRS dgetrs_  
# define  DGETRF dgetrf_  
# define  DGETRI dgetri_  
# define  DGEMM  dgemm_   

int dgetrf_(int *M, int *N, double *A, int *LDA, 
            int *iPiv, int *INFO);

int dgetri_(int *N, double *A, int *LDA, 
            int *iPiv, double *Work, int *WORKL, int *INFO);


int run_lapack(double *matrixWork, int*intWork, double*data)
{
    int n = DIM;
    int ldA = n;
    int info;
    double *Wptr = matrixWork;
    double *Aptr = data;
    int workSize = sizeDoubleWork;
    
    int *iPIV = intWork;
    
    DGETRF(&n,&n,Aptr,&ldA,iPIV,&info);
    if (info != 0) 
      return -abs(info);
    DGETRI(&n,Aptr,&ldA,iPIV,Wptr,&workSize,&info);
}

void print_square(double *A, int n)
{
  for (int i=0; i < n; i++) {
    for (int j=0; j < n; j++)
      printf("  %lf\t", *A++);
    printf("\n");
  }
  printf("\n");
}

int main()
{
  clock_t start, end;
  double cpu_time_used;
  const int n_iter = 10000;

  double W[36];
  int    I[36];
  double A[36] = {
     1.04658065, 1.88413457, 1.47243139, 1.96935254, 1.89419533, 1.58152177,
     1.05004769, 1.31960422, 1.57403921, 1.75055204, 1.5863269 , 1.22141315,
     1.80781111, 1.31407313, 1.17419899, 1.90093874, 1.74648076, 1.48885884,
     1.86851531, 1.86358184, 1.17527866, 1.56880509, 1.44079709, 1.94121331,
     1.58295955, 1.37394298, 1.60141402, 1.57634236, 1.83577456, 1.62298993,
     1.62800519, 1.72560766, 1.96772426, 1.87271493, 1.09485162, 1.27847435,
  };
  start = clock();
  for (int i=0; i<n_iter; i++)
    run_lapack(W, I, A);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Time: %lf\n", cpu_time_used); 
  print_square(A, DIM);

  int ok_flag;
  double B[36] = {
     1.04658065, 1.88413457, 1.47243139, 1.96935254, 1.89419533, 1.58152177,
     1.05004769, 1.31960422, 1.57403921, 1.75055204, 1.5863269 , 1.22141315,
     1.80781111, 1.31407313, 1.17419899, 1.90093874, 1.74648076, 1.48885884,
     1.86851531, 1.86358184, 1.17527866, 1.56880509, 1.44079709, 1.94121331,
     1.58295955, 1.37394298, 1.60141402, 1.57634236, 1.83577456, 1.62298993,
     1.62800519, 1.72560766, 1.96772426, 1.87271493, 1.09485162, 1.27847435
  };

  start = clock();
  for (int i=0; i<n_iter; i++)
    cmx_inv5(B, B, &ok_flag);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Time: %lf\n", cpu_time_used); 
  print_square(B, DIM);

  cmx_inv(B, B, &ok_flag);

#if 0
  double C[36] = {
     1.04658065, 1.88413457, 1.47243139, 1.96935254, 1.89419533, 1.58152177,
     1.05004769, 1.31960422, 1.57403921, 1.75055204, 1.5863269 , 1.22141315,
     1.80781111, 1.31407313, 1.17419899, 1.90093874, 1.74648076, 1.48885884,
     1.86851531, 1.86358184, 1.17527866, 1.56880509, 1.44079709, 1.94121331,
     1.58295955, 1.37394298, 1.60141402, 1.57634236, 1.83577456, 1.62298993,
     1.62800519, 1.72560766, 1.96772426, 1.87271493, 1.09485162, 1.27847435
  };

  start = clock();
  for (int i=0; i<n_iter; i++)
    cmx_inv6_v3(C, C, &ok_flag);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Time: %lf\n", cpu_time_used); 
  print_square(C, DIM);
#endif
}
