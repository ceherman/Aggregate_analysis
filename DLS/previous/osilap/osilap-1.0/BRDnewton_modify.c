#include "ilt.h"
/* ####################################################################### */
/* BRD法で最適化 */
/* 行列サイズは以下のとおり */
/* kz:    M*N (const) */
/* mz:    M*1 (const) */
/* KH:    M*N kz*G  */
/* G:     M*M (diagonal)  */
/* aI:    M*M alpha*(diagonal) */
/* ppchi: M*M */
/* cr:    M*1 */
/* w1:    M*1 */
/* pchi:  M*1 */
/* H:     N*N (diagonal) */
/* h:     N*1 = kz*cr */
/* ####################################################################### */
double BRDnewton_modify(const double *kz, const double *mz, 
		      const double *aI,  double *cr,
		      double *H, HEAD *head, const int *m, const int *n){
  const int incr = 1;
  const int  mm = (*m)*(*m);
  const double zero = 0.0;
  const double one = 1.0;
  const double one_ = -1.0;
  int i, info, step, m_step = 1;
  int *ipiv, lwork = 4*(*m);
  double gamma=1.0, pchi_norm, chi, chi_p=1e300;
  double *p_cr, *p_pchi, *p_ppchi;
  double mz_norm = dnrm2_(m, mz, &incr);
  double *KH, *G, *ppchi;
  double *w1, *dcr, *pchi, *h;

  ipiv= malloc(sizeof(int)*(*n));
  KH = malloc(sizeof(double)*(*m)*(*n));
  G = malloc(sizeof(double)*(*m)*(*m));
  ppchi = malloc(sizeof(double)*(*m)*(*m));
  w1 = malloc(sizeof(double)*(*m));
  dcr = malloc(sizeof(double)*(*m));
  pchi = malloc(sizeof(double)*(*m));
  h = malloc(sizeof(double)*(*n));
  p_cr = malloc(sizeof(double)*(*m));
  p_pchi = malloc(sizeof(double)*(*m));
  p_ppchi = malloc(sizeof(double)*(*m)*(*m));

  double *work, *wi, *wr, *vl, *vr, lambda;
  work = malloc(sizeof(double)*lwork);
  wr = malloc(sizeof(double)*(*m));
  wi = malloc(sizeof(double)*(*m));
  vl = malloc(sizeof(double)*(*m)*(*m));
  vr = malloc(sizeof(double)*(*m)*(*m));

  for (step=0;; step++){
    /* hi = np.dot(k0.T, cr) */
    dgemv_("T", m, n, &one, kz, m, cr, &incr, &zero, h, &incr);
    /* H (diagonal)の計算 */
    for (i=0; i<*n; i++){
      if (h[i]>0) H[i*(*n+1)]=1.0;
      else H[i*(*n+1)]=0.0;
    }
    /* Gr = kz*H*kzT */
    dgemm_("N", "N", m, n, n, &one, kz, m, H, n,  &zero, KH, m, &info);
    dgemm_("N", "T", m, m, n, &one, KH, m, kz, m, &zero, G, m, &info);

    /* chi=0.5*cr*(G+aI)*cr - cr*mz */
    daxpy_(&mm, &one, aI, &incr, G, &incr); /* G = G+aI */
    dcopy_(&mm, G, &incr, ppchi, &incr);    /* G = chi'' */
    dgemv_("N", m, m, &one, G, m, cr, &incr, &zero, w1, &incr);/* w1=(G+aI)cr */
    chi = 0.5*ddot_(m, cr, &incr, w1, &incr);
    chi = chi - ddot_(m, cr, &incr, mz, &incr);
    /* chi' = (G+alphaI)cr-mz */
    dcopy_(m, mz, &incr, pchi, &incr);    /* cr' = mz */
    dscal_(m, &one_, pchi, &incr);        /* cr' = -mz */
    daxpy_(m, &one, w1, &incr, pchi, &incr); /* cr' = w1 -mz */
    pchi_norm = dnrm2_(m, pchi, &incr);

    /* chi'とchi''からdcr=w1=-ppchi/pchiを求める */
    if (chi_p > chi){
      dcopy_(m, cr, &incr, p_cr, &incr);         /* crを保存しておく */
      dcopy_(&mm, ppchi, &incr, p_ppchi, &incr); /* chi''を保存しておく */
      dcopy_(m, pchi, &incr, p_pchi, &incr);     /* chi'を保存しておく */
      dgesv_(m, &incr, ppchi, m, ipiv, pchi, m, &info); /* phiにdcrが入る */
      dcopy_(m, pchi, &incr, dcr, &incr);         /* dcrを保存しておく */
      gamma = -1;
      daxpy_(m, &gamma, dcr, &incr, cr, &incr); /* cr - dcr */
    }
    else { 
      dcopy_(m, p_cr, &incr, cr, &incr);         /* crを元に戻す */
      dcopy_(m, p_pchi, &incr, pchi, &incr);     /* cr'を元に戻す */
      dcopy_(&mm, p_ppchi, &incr, ppchi, &incr); /* cr''を元に戻す */
      /* dcr方向へのステップを減らしてcrを作る */
      gamma = -pow(0.5, (double)m_step);
      daxpy_(m, &gamma, dcr, &incr, cr, &incr); /* cr - dcr */
      m_step++;
      /* 100ステップでchiが減少しなかったら越えたらppchiの正定値を疑う */
      /* Levenberg-Marquardtの方法。固有値を求めて正定値にする */
      if (m_step > 10000){
	printf("Levenberg-Marquardt...\n");
	dgeev_("N","N", m, ppchi, m, wr, wi, vl, m, vr, m, work, &lwork, &info);
	dcopy_(&mm, p_ppchi, &incr, ppchi, &incr); /* cr''を元に戻す */
	lambda = dnrm2_(&mm, ppchi, &incr);
	lambda = wr[idamax_(m, wr, &incr)-1] * lambda;
	if (lambda < 0) lambda = -lambda;
	for (i=0; i<*m; i++) ppchi[i*(*m+1)] = ppchi[i*(*m+1)] + lambda;
	dgesv_(m, &incr, ppchi, m, ipiv, pchi, m, &info); /* phiにdcrが入る */
	dcopy_(m, pchi, &incr, dcr, &incr);
	m_step = 1;
      }
      continue;
    }
    gamma = 1.0;
    m_step = 0;
    chi_p = chi;
    /* 途中経過を表示 */
    printf("%4d: ", step+1);
    /* printf("Dchi: %10.3e", delta_chi); */
    printf(" chi: %10.4g", chi);
    printf(" chi': %10.4g", pchi_norm);
    /* printf(" res : %8.3e", dnrm2_(m, w1, &incr)); */
    printf(" alpha: %10.4e", aI[0]);
    printf(" loop: %d\n", head->alpha_loop+1);
    fflush(stdout);
    /* cr'で収束チェック */
    /* if (dnrm2_(m, w1, &incr) < head->newton_tol) break; */
    if (pchi_norm/mz_norm < head->newton_tol) {
      dcopy_(m, mz, &incr, w1, &incr);
      dgesv_(m, &incr, p_ppchi, m, ipiv, w1, m, &info); /* phiにdcrが入る */
      return dnrm2_(m, w1, &incr);
    }
    if (step == head->newton_loop_max) {
      printf("Iteration reached maximum newton loop.\n");
      break;
    }
  }
  free(work);
  free(ipiv);
  free(KH);
  free(G);
  free(ppchi);
  free(w1);
  free(dcr);
  free(pchi);
  free(h);
  free(p_cr);
  free(p_pchi);
  free(p_ppchi);
  return 0;
}
