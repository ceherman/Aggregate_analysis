#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct{
  double newton_tol;
  double alpha_tol;
  double sfactor;
  double alpha0;
  double sigma0;
  int confirm;
  int alpha_loop_max;
  int newton_loop_max;
  int TXnum, TYnum;
  int DXnum, DYnum;
  int dimension;
  double TXmin, TXmax;
  double TYmin, TYmax;
  int Atype;
  char TXspace[10];
  char TYspace[10];
  char TXtype[10];
  char TYtype[10];
  char infile[1024];
  char outbase[1024];
  int newton_loop;
  int alpha_loop;
  /* s-curve用パラメータ */
  double s_taget;
  double alpha_start;
  double alpha_end;
  int alpha_num;
  double *alpha_list;
  double *sigma_list;

  /* 計算結果用パラメータ */
  double  alpha;
  int s1, s2;
  int N1, N2;
  double sigma;
  FILE* infp;
}HEAD;

typedef struct{
  double *a;
  int m, n;
}MATRIX;

typedef struct{
  MATRIX u; /* tau1のベクトル */
  MATRIX v; /* tau2のベクトル */
  MATRIX m; /* 観測データの行列 u.m*v.m */
  double mmax;  /* 観測データの最大値 */
}INDATA;

double ddot_();
double dnrm2_();
void dcopy_();
void dgemm_();
void dgesvd_();
void dgesv_();
void dgeev_();
void dsysv_();
void dgemv_();
void nnls_c();
void daxpy_();
void dscal_();
int idamax_();
int SetArgment();
void Transpose();
void OutputDistr();
void OutputDecay();
void OutputDistr1D();
void OutputDecay1D();
void Error();
void Confirm();
void MakeKernel();
void MultiplyMatrix();
void Lexicographic();
void ReadData();
void ReadData1D();
void SetDataSize();
double BRDnewton_modify();
void BRDcg();
void BRDAlphaLoop();
void ScurveAlphaLoop();
void CalcSim();
double CalcSigma();
void OutputScurve();
