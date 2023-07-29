#include "ilt.h"

/* ####################################################################### */
/* 最適化したfから圧縮無しのσ^2を求める */
/* ####################################################################### */
double CalcSigma(MATRIX *K1, MATRIX *K2, double *f, HEAD *head, INDATA *d){
  int i, incr = 1;
  int m1 = K1->m, m2 = K2->m;
  int n1 = K1->n, n2 = K2->n;
  int nn = n1*n2, mm = m1*m2;
  int info;
  double one=1.0, one_=-1.0, zero = 0.0, chi, sigma=0;
  double *KF, *F, *M;
  div_t temp;

  if (head->dimension == 2){
    KF = malloc(sizeof(double)*m1*n2);
    F = malloc(sizeof(double)*nn);
    M = malloc(sizeof(double)*m1*m2);
    /* fから行列Fを生成 */
    for (i=0; i<nn; i++){
      temp = div(i, n2);
      F[temp.rem*n1 + temp.quot] = f[i];
    }

    /* K1*Fの計算 */
    dgemm_("N", "N", &m1, &n2, &n1, &one, K1->a, &m1, 
	   F, &n1, &zero, KF, &m1, &info);
    /* (K1*F)*K2Tの計算 */
    dgemm_("N", "T", &m1, &m2, &n2, &one, KF, &m1, 
	   K2->a, &m2, &zero, M, &m1, &info);
    /* M = -Mexp + Mcalの計算 */  
    daxpy_(&mm, &one_, d->m.a, &incr, M, &incr);
    /* ||Mcal - Msim|| */
    chi = dnrm2_(&mm, M, &incr);
    /* sigma = σ^2 = ||KF-M||^2/m */
    sigma = chi*chi/mm;
    free(KF);
    free(F);
    free(M);
  }
  if (head->dimension == 1){
    M = malloc(sizeof(double)*m1);
    dgemv_("N", &m1, &n1, &one, K1->a, &m1, f, &incr, &zero, M, &incr); 
    /* M = -Mexp + Mcalの計算 */  
    daxpy_(&m1, &one_, d->m.a, &incr, M, &incr);
    chi = dnrm2_(&m1, M, &incr);
    /* sigma = σ^2 = ||KF-M||^2/m */
    sigma = chi*chi/m1;
    free(M);
  }
  return sigma;
}

/* ####################################################################### */
/* Scurve法でalpha_optを求める */
/* ####################################################################### */
void ScurveAlphaLoop(MATRIX *K0, MATRIX *Mz, HEAD *head, double *f,
		     MATRIX *k1, MATRIX *k2, INDATA *d){
  int m = K0->m;
  int n = K0->n;
  int i, j, pnt = 0, incr = 1;
  int  alpha_num=head->alpha_num;
  double *kz, *mz, *cr, *aI, *H;
  double alpha, sigma;
  double zero = 0.0, one = 1.0;
  double alpha_s = head->alpha_start, alpha_e = head->alpha_end;
  double s_target=head->s_taget, alpha_log;
  double *alpha_list, *sigma_list, *grad;
  double *t = NULL;
  
  kz = K0->a;
  mz = Mz->a;
  cr = malloc(sizeof(double)*m);
  aI = malloc(sizeof(double)*m*m);
  H = malloc(sizeof(double)*n*n);
  alpha_list = (double *)malloc(sizeof(double)*alpha_num);
  sigma_list = (double *)malloc(sizeof(double)*alpha_num);
  grad = malloc(sizeof(double)*(alpha_num-1));
  /* alphaのグリッド生成 */
  for (i=0; i<alpha_num; i++){
    alpha_log = alpha_s + (double)i*(alpha_e-alpha_s)/(alpha_num-1);
    alpha_list[i] = pow(10, alpha_log);
  }
  
  /* alphaI (diagonal) とH (diagonal)の初期化 */
  for (i=0; i<(n*n); i++) H[i] = 0.0;
  for (i=0; i<m*m; i++) aI[i] = 0.0;

  /* ここからalphaのloop */
  for (j=alpha_num-1; j>=0; j--){
    /* alphaI (diagonal)の初期化*/
    for (i=0; i<m; i++) aI[i*(m+1)] = alpha_list[j];
    /* BRDでcrの最適化 */
    BRDnewton_modify(kz, mz,  aI, cr,  H, head, &m, &n);
    /* crからfの計算 f= kz'*cr */
    dgemv_("T", &m, &n, &one, kz, &m, cr, &incr, &zero, f, &incr);
    /* fi>0 だけが解 */
    for (i=0; i<n; i++) (f[i] < 0) ? (f[i] = 0.0) : (f[i] = f[i]);
    /* 圧縮前のカーネルを使って圧縮ノイズの標準偏差を求める */
    sigma_list[j] = CalcSigma(k1, k2, f, head, d); 
    head->alpha_loop++;
  }
  
  for (i=0; i<alpha_num-1; i++){
    grad[i] = (log10(sigma_list[i+1]) - log10(sigma_list[i]))/
              (log10(alpha_list[i+1]) - log10(alpha_list[i]));
    if (grad[i] > s_target && pnt == 0) pnt = i;
  }
  double alpha_0, alpha_1, grad_0, grad_1;
  double sigma_0, sigma_1;
  int step;
  alpha_0 = alpha_list[pnt-1];
  sigma_0 = sigma_list[pnt-1];
  alpha_1 = alpha_list[pnt+1];
  sigma_1 = sigma_list[pnt+1];
  for (step=1;; step++){
    pnt = alpha_num+step;
    alpha = alpha_0 + 0.5*(alpha_1 - alpha_0);
    if ((t = (double *)realloc(alpha_list, pnt*sizeof(double))) != NULL){
      alpha_list = t;
      alpha_list[alpha_num+step-1] = alpha;
    }
    for (i=0; i<m; i++) aI[i*(m+1)] = alpha;
    BRDnewton_modify(kz, mz,  aI, cr,  H, head, &m, &n);
    dgemv_("T", &m, &n, &one, kz, &m, cr, &incr, &zero, f, &incr);
    for (i=0; i<n; i++) (f[i] < 0) ? (f[i] = 0.0) : (f[i] = f[i]);
    sigma = CalcSigma(k1, k2, f, head, d); 
    if ((t = (double *)realloc(sigma_list, pnt*sizeof(double))) != NULL){
      sigma_list = t;
      sigma_list[pnt-1] = sigma;
    }
    grad_0 = (log10(sigma) - log10(sigma_0))/(log10(alpha) - log10(alpha_0));
    grad_1 = (log10(sigma_1) - log10(sigma))/(log10(alpha_1) - log10(alpha));
    if (fabs(s_target - grad_0) <  fabs(s_target - grad_1) ){
      alpha_1 = alpha;
      sigma_1 = sigma;
    }
    else {
      alpha_0 = alpha;
      sigma_0 = sigma;
    }
    head->alpha = alpha;
    if ((alpha_1 - alpha_0) < head->alpha_tol) break;
    printf("delta alpha = %g\n", alpha_1 - alpha_0);
    printf("grad0 = %g\n", grad_0);
    printf("grad1 = %g\n", grad_1);
    head->alpha_loop++;
  }

  head->alpha_list = malloc(pnt*sizeof(double));
  head->sigma_list = malloc(pnt*sizeof(double));
  head->alpha_num = pnt;
  for (i=0; i<pnt; i++){
    head->alpha_list[i] = alpha_list[i];
    head->sigma_list[i] = sigma_list[i];
  }

  free(cr);
  free(aI);
  free(H);
  free(alpha_list);
  free(sigma_list);
  free(grad);
}

void BRDAlphaLoop(MATRIX *K0, MATRIX *Mz, HEAD *head, double *f,
		  MATRIX *k1, MATRIX *k2, INDATA *d){
  int m = K0->m;
  int n = K0->n;
  int i, incr = 1;
  double *kz, *mz;
  double *cr, *cr_n;
  double *aI, *H;
  double zero = 0.0, one = 1.0, one_ = -1.0;
  double alpha_opt, cr_norm, sigma, sigma_z;
  
  if (head->alpha0  == 0) {
    head->alpha = 0.0;
    return;
  }

  kz = K0->a;
  mz = Mz->a;
  cr = malloc(sizeof(double)*m);
  cr_n = malloc(sizeof(double)*m);
  aI = malloc(sizeof(double)*m*m);
  H = malloc(sizeof(double)*n*n);
  
  /* 初期のfからcrを求める cr=-(Kf-m)/alpha */
  dgemv_("N", &m, &n, &one, kz, &m, f, &incr, &zero, cr, &incr);
  daxpy_(&m, &one_, mz, &incr, cr, &incr);
  alpha_opt = -1.0/head->alpha0; 
  dscal_(&m, &alpha_opt, cr, &incr);
  head->alpha = head->alpha0;

  /* alphaI (diagonal) とH (diagonal)の初期化 */
  for (i=0; i<(n*n); i++) H[i] = 0.0;
  for (i=0; i<m*m; i++) aI[i] = 0.0;
  /* for (i=0; i<m; i++) cr[i] = 0.0; */
  /* ここからalphaのloop */
  for (head->alpha_loop=0;; head->alpha_loop++){
    /* alphaI (diagonal)の初期化*/
    for (i=0; i<m; i++) aI[i*(m+1)] = head->alpha;
    /* BRDでcrの最適化 */
    cr_norm = BRDnewton_modify(kz, mz,  aI, cr,  H, head, &m, &n);
    /* crからfの計算f f= kz'*cr,  fi>0 だけが解 */
    dgemv_("T", &m, &n, &one, kz, &m, cr, &incr, &zero, f, &incr);
    for (i=0; i<n; i++) (f[i] < 0) ? (f[i] = 0.0) : (f[i] = f[i]);
    /* 圧縮前のカーネルを使って圧縮ノイズの標準偏差を求める */
    sigma = sqrt(CalcSigma(k1, k2, f, head, d));
    printf("uncompressed sigma = %g\n", sigma);
    /* 圧縮のカーネルを使ってcrを求める */    
    dgemv_("N", &m, &n, &one, kz, &m, f, &incr, &zero, cr_n, &incr);
    daxpy_(&m, &one_, mz, &incr, cr_n, &incr); /* kz*f - mz */
    sigma_z = dnrm2_(&m, cr_n, &incr)/sqrt(m);
    printf("compressed sigma = %g\n", sigma_z);

    if (head->sigma0 != 0.0) sigma = head->sigma0;
    printf("c_norm = %g\n", cr_norm);

    switch(head->Atype){
    case 0:
      alpha_opt = sigma/sigma_z;
      break;
    case 1:
      alpha_opt = sigma_z/dnrm2_(&n, f, &incr)/(sqrt(1-sqrt(2*m)/m));
      break;
    case 2:
      alpha_opt = sigma_z/dnrm2_(&n, f, &incr);
      break;
    case 3:
      alpha_opt = sigma/dnrm2_(&n, f, &incr)/(sqrt(1-sqrt(2*m)/m));
      break;
    }

    /* alphaの収束チェック */
    if (fabs(alpha_opt - head->alpha)/head->alpha < head->alpha_tol)
      break;
    /* alphaをmixingしてupdate */
    if (head->alpha_loop > 1) 
      head->alpha = alpha_opt;
    else 
    head->alpha_loop++;
    if (head->alpha_loop == head->alpha_loop_max) {
      printf("Iteration reached maximum alpha loop.\n");
      break;
    }
  }
  free(cr_n);
  free(cr);
  free(aI);
  free(H);
}
