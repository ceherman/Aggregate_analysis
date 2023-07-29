#include "ilt.h"

/* ####################################################################### */
/* Setting Parameters */
/* default values */
/* ####################################################################### */
void SetParameter(HEAD *head){
  head->newton_tol = 1e-8;
  head->alpha_tol = 1e-4;
  head->confirm = 0;
  head->dimension = 2;
  strcpy(head->outbase, "");
  head->Atype = 0;
  head->sfactor = 1.0e-4;
  head->alpha0 = 0.1;
  head->alpha_loop_max = 5000;
  head->newton_loop_max = 5000;
  head->sigma0 = 0;
  head->DXnum = 0;
  head->DYnum = 0;
  head->TXnum = 50; 
  head->TXmin = 1.0e-4;
  head->TXmax = 10;
  strcpy(head->TXspace, "log");
  strcpy(head->TXtype, "T1");

  head->TYnum = 50;
  head->TYmin = 1.0e-4;
  head->TYmax = 10;
  strcpy(head->TYspace, "log");
  strcpy(head->TYtype, "T2");

  head->newton_loop = 0;
  head->alpha_loop  = 0;

  head->s_taget = 0.0;
  head->alpha_start = -4; /* log10 space */
  head->alpha_end = 6;    /* log10 space */
  head->alpha_num = 20;
}

/* ####################################################################### */
/* Transpose Matrix */
/* ####################################################################### */
void Transpose(MATRIX *A, MATRIX *B){
  int x;
  int i, j, p;
  div_t temp;

  B->m = A->n;
  B->n = A->m;
  B->a = (double*)malloc(sizeof(double)*B->m*B->n);
  for (x=0; x<A->m*A->n; x++){
    temp = div(x, A->m);
    i = temp.rem;  j = temp.quot;
    p = i*B->m + j;
    B->a[p] = A->a[x];
  }
}

/* ####################################################################### */
/* Kronecker (Tensor) product */
/* ####################################################################### */
void Kronecker(MATRIX *k1, MATRIX *k2, MATRIX *k0){
  int i, j, k, l;
  int row, col;
  int  p;
  int x, y;
  div_t temp;
  
  k0->m = k1->m*k2->m;
  k0->n = k1->n*k2->n;
  k0->a = (double*)malloc(sizeof(double)*k0->m*k0->n);

  for (x=0; x<k1->m*k1->n; x++){
    temp = div(x, k1->m);
    i = temp.rem;  j = temp.quot;
    for (y=0; y<k2->m*k2->n; y++){
      temp = div(y, k2->m);
      k = temp.rem;  l = temp.quot;
      row = i*k2->m + k;
      col = j*k2->n + l;
      p = col * k0->m + row;
      k0->a[p] = k1->a[x] * k2->a[y];
    }
  }
}

/* ####################################################################### */
/* Lexicographic row major */
/* ####################################################################### */
void Lexicographic(MATRIX *A, MATRIX *Al){
  int i, p;
  const int incr = 1;
  const int mn = A->m*A->n;
  div_t temp;

  Al->a = (double*)malloc(sizeof(double)*A->m*A->n);
  Al->m = A->m * A->n;
  Al->n = 1;

  dcopy_(&mn, A->a, &incr, Al->a, &incr);
  for (i=0; i<mn; i++){
    temp = div((int)i, (int)A->n);
    p = temp.rem * A->m + temp.quot;
    Al->a[i] = A->a[p];
  }
}

/* ####################################################################### */
/* 行列の掛け算 C = A.a*B.a + beta*C.a */
/* alpha=beta=1 */
/* ####################################################################### */
void MultiplyMatrix(MATRIX *A, MATRIX *B, MATRIX *C, 
		    char aopt[], char bopt[]){
  int m, k, n; 
  int lda, ldb, ldc;
  double alpha = 1.0;
  double beta  = 0.0;
  int  info;  

  if (strcmp(aopt, "N")==0){
    m = A->m;   k = A->n;
    lda = A->m; ldb = B->m; ldc = A->m;
    C->m = A->m;
  }
  if (strcmp(aopt, "T")==0){
    m = A->n;   k = A->m;
    lda = A->m; ldb = B->m; ldc = A->n;
    C->m = A->n;
  }
  if (strcmp(bopt, "N")==0){
    n = B->n; ldb = B->m;
    C->n = B->n;
  }
  if (strcmp(bopt, "T")==0){
    n = B->m; ldb = B->m;
    C->n = B->m;
  }

  C->a = (double*)malloc(sizeof(double)*C->m*C->n);
  dgemm_(aopt, bopt, &m, &n, &k, &alpha, A->a, &lda, 
	 B->a, &ldb, &beta, C->a, &ldc, &info);
}

/* ####################################################################### */
/* 配列の出力 */
/* ####################################################################### */
void PrintList(const double *x, const int m){
  int i;

  printf("Vector size: %d\n", m);
  for (i=0; i<m; i++){
    printf("%16.8e\n", x[i]);
  }
  exit(0);
}

/* ####################################################################### */
/* 行列の出力 */
/* ####################################################################### */
void PrintMatrix(const double *mat, const int m, const int n, 
		 const int opt){
  int i, j;
  
  printf("Matrix size: %d %d\n", m, n);
  for (i=0; i<m; i++){
    if ( 2 < i && i < m-3 && opt==1) continue;
    for (j=0; j<n; j++){
      if ( 2 < j && j < n-3 && opt==1) continue;
      printf("%16.8e", mat[j+n*(i)]);
      if (opt==0) printf("\n");
    }
    if (opt==1) printf("\n");
  }
  exit(0);
}

/* ####################################################################### */
/* 緩和分布の離散化 */
/* ####################################################################### */
void MakeKspace(double start, double end, int points, 
		MATRIX *x, char *type){
  int i;
  double dx;

  x->m = points;
  x->n = 1;
  x->a = (double*)malloc(sizeof(double)*points);
  if (strcmp(type, "log")==0){
    start = log10(start);
    end = log10(end);
    dx = (end - start)/(double)(points-1);
    for (i=0; i<points; i++) x->a[i] =  pow(10, start+i*dx);
  }
  else if (strcmp(type, "linear")==0){
    start = start;
    end = end;
    dx = (end - start)/(double)(points-1);
    for (i=0; i<points; i++) x->a[i] =  start+i*dx;
  }
  else {
    printf("Kernel space %s is strange.", type);
    Error();
  }
}

/* ####################################################################### */
/* kernel生成 */
/* ####################################################################### */
void MakeKernel(MATRIX *t, MATRIX *x, char *type, MATRIX *k){
  int i, j, c=0;

  if ((strcmp(type, "T1inf") == 0) || 
      (strcmp(type, "T2inf") == 0) ||
      (strcmp(type, "T2ginf") == 0) ){
    k->a = (double*)malloc(sizeof(double)*t->m*(x->m+1));
    k->m = t->m;
    k->n = x->m + 1;
  }
  else {
    k->a = (double*)malloc(sizeof(double)*t->m*x->m);
    k->m = t->m;
    k->n = x->m;
  }
  if (strcmp(type, "T1") == 0 || 
      strcmp(type, "T1inf") == 0){
    for (i=0; i<x->m; i++){
      for (j=0; j<t->m; j++){
	k->a[c] = 1.0-2.0*exp(-t->a[j]/x->a[i]);
	c++;
      }
    }
  }
  if (strcmp(type, "T2") == 0 || 
      strcmp(type, "T2inf") == 0 || 
      strcmp(type, "T1p") == 0){
    for (i=0; i<x->m; i++){
      for (j=0; j<t->m; j++){
	k->a[c] = exp(-t->a[j]/x->a[i]);
	c++;
      }
    }
  }
  if (strcmp(type, "D") == 0 ){
    for (i=0; i<x->m; i++){
      for (j=0; j<t->m; j++){
	k->a[c] = exp(-t->a[j]*x->a[i]);
	c++;
      }
    }
  }
  if (strcmp(type, "T2g") == 0 ||
      strcmp(type, "T2ginf") == 0){
    for (i=0; i<x->m; i++){
      for (j=0; j<t->m; j++){
	k->a[c] = exp(-0.5*(t->a[j]/x->a[i])*(t->a[j]/x->a[i]));
	c++;
      }
    }
  }
  
  /* base lineの継ぎ足し */
  if (strcmp(type, "T1inf") == 0 || 
      strcmp(type, "T2inf") == 0||
      strcmp(type, "T2ginf") == 0){
    for (j=0; j<t->m; j++){
      k->a[c] = 1.0;
      c++;
    }
  }
}

/* ####################################################################### */
/* SVD compress */
/* A = U * S * Vt */
/* ####################################################################### */
void SvdCompress(MATRIX *A, MATRIX *uz, MATRIX *sz, MATRIX *vz, 
		 double sfactor){
  int i, j;
  const int incr = 1;
  const int  m = A->m, n= A->n, mn = A->m*A->n;
  const int lwork=5*m*n;
  int info;
  double *U, *Vt, *Aa, *work;
  MATRIX S;

  (m < n) ? (S.m = A->m) : (S.m = A->n);
  S.n = 1;
  Aa = malloc(sizeof(double)*m*n);
  U = malloc(sizeof(double)*m*m);
  S.a = malloc(sizeof(double)*m*n);
  Vt = malloc(sizeof(double)*n*n);
  work = malloc(sizeof(double)*lwork);
  dcopy_(&mn, A->a, &incr, Aa, &incr);
  dgesvd_("A", "A", &m, &n, Aa, &m, S.a, 
	  U, &m, Vt, &n, work, &lwork, &info);
  if (info != 0){
    printf("svd is fault!!! info=%d\n", info);
    exit(1);
  }

  /* 切り捨てずに採用するランク数の計算 */
  for (i=0; S.a[i]>S.a[0]*sfactor; i++);
  sz->m = i; sz->n = i;
  sz->a = (double*)malloc(sizeof(double) * sz->m * sz->n);
  uz->m = m; uz->n = i;
  uz->a = (double*)malloc(sizeof(double) * uz->m * uz->n);
  vz->m = i; vz->n = n;
  vz->a = (double*)malloc(sizeof(double) * vz->m * vz->n);

  /* 圧縮した行列の生成 */  
  for (i=0; i < (sz->m*sz->n); i++) sz->a[i] = 0.0;      /* 初期化 */
  for (i=0; i < sz->n; i++) sz->a[i*(sz->n+1)] = S.a[i]; /* diagonalだけ */
  for (i=0; i < (uz->n*uz->m); i++) uz->a[i] = U[i];
  for (i=0; i < (vz->n); i++){
    for (j=0; j < vz->m; j++){
      vz->a[j+i*vz->m] = Vt[j+i*n]; 
    }
  }
  free(S.a);
  free(Aa);
  free(U);
  free(Vt);
  free(work);
}

/* ####################################################################### */
/* NNLSで初期値を決める */
/* ####################################################################### */
void Calcf0(MATRIX *k0, MATRIX *mz, double *f0){
  double rnorm;
  const int incr = 1;
  const int m=k0->m;
  const int n=k0->n;
  const int mn=k0->m*k0->n;
  int index[k0->n];
  int mode;
  double *mn_mem, *m_mem, *m_mem2, *n_mem;

  mn_mem = malloc(sizeof(double)*k0->m*k0->n);
  m_mem = malloc(sizeof(double)*k0->m);
  m_mem2 = malloc(sizeof(double)*k0->m);
  n_mem = malloc(sizeof(double)*k0->n);

  dcopy_(&mn, k0->a, &incr, mn_mem, &incr);
  dcopy_(&m, mz->a, &incr, m_mem, &incr);
  nnls_c(mn_mem, &m, &m, &n, m_mem, f0, &rnorm, 
	 n_mem, m_mem2, index, &mode);
  free(mn_mem);  
  free(m_mem);
  free(m_mem2);   
  free(n_mem);
}


void CalcSim1D(MATRIX *k, double *f, double *baseline, MATRIX *Msim,
	       INDATA *d, HEAD *head){
  int incr = 1;
  double alpha = 1.0, beta = 0.0;
  double *deltaM;
  double one_ = -1.0;

  Msim->a = malloc(sizeof(double)*k->m);
  deltaM = malloc(sizeof(double)*k->m);
  dgemv_("N", &k->m, &k->n, &alpha, k->a, &k->m, f, 
	 &incr, &beta, Msim->a, &incr);
  if (strcmp(head->TXtype, "T2inf")==0){
    baseline[0] = f[k->n-1];
  }

  dcopy_(&k->m, Msim->a, &incr, deltaM, &incr);
  daxpy_(&k->m, &one_, d->m.a, &incr, deltaM, &incr);
  head->sigma = dnrm2_(&k->m, deltaM, &incr)/sqrt(k->m);
  free(deltaM);

}

void CalcSim(MATRIX *k1, MATRIX *k2, double *f, double *baseline, MATRIX *Msim, 
	     INDATA *d, HEAD *head){
  int i, j, k=0, mn = k1->m*k2->n;
  int mm = k1->n * k2->n;
  MATRIX F, KF;
  div_t temp;
  double *deltaM;
  double one_ = -1.0;
  int incr = 1;

  F.m = k1->n;
  F.n = k2->n;
  F.a = malloc(sizeof(double)*F.m*F.n);
  deltaM = malloc(sizeof(double)*mm);
  for (i=0; i<F.n; i++) baseline[i] = 0.0;
  mn = F.m * F.n;
  for (i=0; i<mn; i++){
    temp = div(i, F.n);
    F.a[temp.rem*F.m + temp.quot] = f[i];
  }
  if (strcmp(head->TXtype, "T1inf") == 0){
    j = 0;
    for (i=0; i<mn; i++){
      temp = div(i, F.n);
      if (temp.quot != (F.m-1)){
	f[j] = F.a[temp.rem*F.m + temp.quot];
	j++;
      }
      else{
	baseline[k] = F.a[temp.rem*F.m + temp.quot];
	k++;
      }
    }
  }
  MultiplyMatrix(k1, &F, &KF, "N", "N");
  MultiplyMatrix(&KF, k2, Msim,  "N", "T");
  dcopy_(&mm, Msim->a, &incr, deltaM, &incr);
  daxpy_(&mm, &one_, d->m.a, &incr, deltaM, &incr);
  head->sigma = dnrm2_(&mm, deltaM, &incr)/sqrt(mm);
  free(F.a);
  free(deltaM);
}

void Copyright(){
  int i;
  printf("\n");
  for (i=0; i<70; i++) printf("-");
  printf("\n");
  printf("1D- and 2D-ILT program, MAY 2014\n");
  printf("Author: T. Ohkubo (Chiba University)\n\n");
  printf("The program can not be distributed without author's permission.\n");
  printf("Use of this program is limited to non-profit or evaluation purposes.");

  printf("\n");
  for (i=0; i<70; i++) printf("-");
  printf("\n");
  printf("\n");

}

int main(int argc, char** argv){
  INDATA d;
  MATRIX k0, k1, k2, s1z, s2z, u1z, u2z, v1z, v2z;
  MATRIX k1z, k2z, mz_tmp1, mz_tmp2, mz, msim;
  MATRIX x, y;      /* T1とT2の離散化データ */
  double *f, *baseline;
  HEAD head;
  Copyright();

  SetParameter(&head);
  SetArgment(argc, argv, &head);
  MakeKspace(head.TXmin, head.TXmax, head.TXnum, &x, head.TXspace);
  MakeKspace(head.TYmin, head.TYmax, head.TYnum, &y, head.TYspace);
  if (head.dimension == 2){
    SetDataSize(&head, &d);
    ReadData(&head, &d);
    Confirm(&head, &d);
    printf("Making kernel...\n");
    MakeKernel(&d.u, &x, head.TXtype, &k1);
    MakeKernel(&d.v, &y, head.TYtype, &k2);
    printf("SVD compressing...\n");
    printf("   Kernel size (K1): %d*%d\n", k1.m, k1.n);
    printf("   Kernel size (K2): %d*%d\n", k2.m, k2.n);
    SvdCompress(&k1, &u1z, &s1z, &v1z, head.sfactor);
    SvdCompress(&k2, &u2z, &s2z, &v2z, head.sfactor);
    printf("Making compress kernel...\n");
    MultiplyMatrix(&s1z, &v1z, &k1z, "N", "N");
    MultiplyMatrix(&s2z, &v2z, &k2z, "N", "N");
    printf("   Kernel size (K1): %d*%d\n", k1z.m, k1z.n);
    printf("   Kernel size (K2): %d*%d\n", k2z.m, k2z.n);
    head.s1 = k1z.m; head.N1 = k1.m;
    head.s2 = k2z.m; head.N2 = k2.m;
    MultiplyMatrix(&u1z, &d.m, &mz_tmp1, "T", "N");
    MultiplyMatrix(&mz_tmp1, &u2z, &mz_tmp2, "N", "N");
    Lexicographic(&mz_tmp2, &mz);
    Kronecker(&k1z, &k2z, &k0);
    /* mz, k0以外のmallocしたメモリをrelease */
    free(s1z.a); free(s2z.a);
    free(v1z.a); free(v2z.a);
    free(u1z.a); free(u2z.a);
    free(k1z.a); free(k2z.a);
    free(mz_tmp1.a); free(mz_tmp2.a);
  }
  if (head.dimension == 1){
    ReadData1D(&head, &d);
    Confirm(&head, &d);
    printf("Making kernel...\n");
    MakeKernel(&d.u, &x, head.TXtype, &k1);
    printf("Making compress kernel...\n");
    SvdCompress(&k1, &u1z, &s1z, &v1z, head.sfactor);
    MultiplyMatrix(&s1z, &v1z, &k0, "N", "N");
    head.s1 = k0.m;     head.N1 = k1.m;
    k2.n = 1;
    printf("   Kernel size (K1): %d*%d\n", k0.m, k0.n);
    MultiplyMatrix(&u1z, &d.m, &mz, "T", "N");
  }
  printf("Calculating initial F using NNLS... \n"); 
  f = malloc(sizeof(double)*k0.n);
  baseline = malloc(sizeof(double)*k2.n);

  Calcf0(&k0, &mz, f);
  
  if (head.s_taget == 0) {
    printf("Start BRD routine...\n");
    BRDAlphaLoop(&k0, &mz, &head, f, &k1, &k2, &d);
  }
  else  {
    printf("Start S-curve routine...\n");
    ScurveAlphaLoop(&k0, &mz, &head, f, &k1, &k2, &d);
  }
  printf("Output results...\n");
  if (head.dimension == 2){
    CalcSim(&k1, &k2, f, baseline, &msim, &d, &head);
    OutputDistr(&head, &x, &y, f, baseline);
    OutputDecay(&head, &d, &x, &y, &msim);
  }
  if (head.dimension == 1){
    CalcSim1D(&k1, f, baseline, &msim, &d, &head);
    OutputDistr1D(&head, &x, f, baseline);
    OutputDecay1D(&head, &d, &x, &msim);
  }
  Copyright();
  return 0;
}
