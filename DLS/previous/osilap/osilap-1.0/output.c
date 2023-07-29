#include "ilt.h"

/* ####################################################################### */
/* Make header information */
/* ####################################################################### */
void OutputHeader(FILE *fp, HEAD *head){
  fprintf(fp, "# -------------- solving condition -------------- #\n");
  fprintf(fp, "#  sfactor:         %e\n", head->sfactor);
  fprintf(fp, "#  initial alpha:   %g\n", head->alpha0);
  fprintf(fp, "#  alpha max loop:  %d\n", head->alpha_loop_max);
  fprintf(fp, "#  alpha tol:       %g\n", head->alpha_tol);
  fprintf(fp, "#  Newton max loop: %d\n", head->newton_loop_max);
  fprintf(fp, "#  Newton tol:      %e\n", head->newton_tol);
  fprintf(fp, "#  A_opt method:    %d\n", head->Atype);
  fprintf(fp, "# --------------- TX information ---------------- #\n");
  fprintf(fp, "#  Data space:      %s\n", head->TXspace);
  fprintf(fp, "#  Data type:       %s\n", head->TXtype);
  fprintf(fp, "#  (TXmin, TXmax, TXpoints)=");
  fprintf(fp, "(%g, %g, %d)\n", head->TXmin, head->TXmax, head->TXnum);
  if (head->dimension == 2){
    fprintf(fp, "# ---------------- TY information --------------- #\n");
    fprintf(fp, "#  Data space:      %s\n", head->TYspace);
    fprintf(fp, "#  Data type:       %s\n", head->TYtype);
    fprintf(fp, "#  (TYmin, TYmax, TYpoints)=");
    fprintf(fp, "(%g, %g, %d)\n", head->TYmin, head->TYmax, head->TYnum);
  }
  fprintf(fp, "# ------------------ input data ----------------- #\n");
  fprintf(fp, "#  Input file:      %s\n", head->infile);
  fprintf(fp, "#  Output basename: %s\n", head->outbase);
  fprintf(fp, "# ---------------- Compress data ---------------- #\n");
  fprintf(fp, "#  K1 original :    %d\n", head->N1);
  fprintf(fp, "#  K1 compress :    %d\n", head->s1);
  if (head->dimension == 2){
    fprintf(fp, "#  K2 original :    %d\n", head->N2);
    fprintf(fp, "#  K2 compress :    %d\n", head->s2);
  }
  fprintf(fp, "# ------------- Calculation results ------------- #\n");
  fprintf(fp, "#  Alpha :          %e\n", head->alpha);
  fprintf(fp, "#  Sigma :          %e\n", head->sigma);
  fprintf(fp, "# ----------------------------------------- #\n");
  fprintf(fp, "\n");
}

/* ####################################################################### */
/* output for distribution data*/
/* ####################################################################### */
void OutputDistr(HEAD *head, MATRIX *x, MATRIX *y, double *f, double *baseline){
  int i, j, k=0;
  char outfile[2048];
  FILE *fp;
  double dmax;
  int incr = 1;
  int mn = x->m * y->m;

  printf("  Distribution data: %s.distr\n", head->outbase);
  if ((strcmp(head->outbase, "-")==0)) fp = stdout;
  else { 
    sprintf(outfile, "%s.distr", head->outbase);
    fp = fopen(outfile, "w");
  }
  OutputHeader(fp, head);
  dmax = fabs(f[idamax_(&mn, f, &incr)-1]);
  fprintf(fp, "# Distribution normalized coeff.: %g\n", dmax);
  fprintf(fp, "# %13s %15s %15s", "TX", "TY", "I\n");
  for (i=0; i<x->m; i++){
    for (j=0; j<y->m; j++){
      fprintf(fp, "%15.8e %15.8e %15.8e\n", x->a[i], y->a[j], f[k]/dmax);
      k++;
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n\n\n");
  fprintf(fp, "# Baseline list\n");
  for (i=0; i<head->TYnum; i++){
    fprintf(fp, "# %15.8e\n", baseline[i]/dmax);    
  }
  fclose(fp);
}

/* ####################################################################### */
/* output for decay data */
/* ####################################################################### */
void OutputDecay(HEAD *head, INDATA *d, 
		 MATRIX *x, MATRIX *y, MATRIX *m){
  int i, j, p=0;
  double t1, t2;
  char outfile[2048];
  MATRIX Iexp, Isim;
  FILE *fp;

  printf("  Decay data: %s%s\n", head->outbase, ".decay");
  Transpose(&d->m, &Iexp);
  Transpose(m, &Isim);
  if ((strcmp(head->outbase, "-")==0)) fp = stdout;
  else { 
    sprintf(outfile, "%s.decay", head->outbase);
    fp = fopen(outfile, "w");
  }
  OutputHeader(fp, head);
  fprintf(fp, "# Signal normalized coeff.: %g\n", d->mmax);
  fprintf(fp, "# %13s %15s %15s %15s %15s\n", 
	  "t1", "t2", "Iexp", "Isim", "Iexp - Isim");
  for (i=0; i<d->u.m; i++){
    for (j=0; j<d->v.m; j++){
      t1 = d->u.a[i];
      t2 = d->v.a[j];
      fprintf(fp, "%15.8e %15.8e %15.8e %15.8e %15.8e\n",
      	      t1, t2, Iexp.a[p], Isim.a[p], Iexp.a[p] - Isim.a[p]);
      p++;
    }
    fprintf(fp, "\n\n");
  }
  fclose(fp);
  if (head->s_taget != 0) OutputScurve(head);
}

/* ####################################################################### */
/* output for distribution data */
/* ####################################################################### */
void OutputDistr1D(HEAD *head, MATRIX *x, double *f, double *baseline){
  int i;
  char outfile[2048];
  FILE *fp;
  double dmax;
  int incr = 1;
  int mn = x->m;

  printf("  Distribution data: %s.distr\n", head->outbase);
  if ((strcmp(head->outbase, "-")==0)) fp = stdout;
  else { 
    sprintf(outfile, "%s%s", head->outbase, ".distr");
    fp = fopen(outfile, "w");
  }
  OutputHeader(fp, head);
  dmax = fabs(f[idamax_(&mn, f, &incr)-1]);
  fprintf(fp, "# Distribution normalized coeff.: %g\n", dmax);
  fprintf(fp, "# %13s %15s \n", "TX", "I");
  for (i=0; i<x->m; i++){
    fprintf(fp, "%15.8e %15.8e\n", x->a[i], f[i]/dmax);
  }
  fprintf(fp, "\n\n\n");
  fprintf(fp, "# Baseline list\n");
  fprintf(fp, "# %15.8e\n", baseline[0]/dmax);    
  if ((strcmp(head->outbase, "-")!=0))   fclose(fp);
}

/* ####################################################################### */
/* output for decay data */
/* ####################################################################### */
void OutputDecay1D(HEAD *head, INDATA *d,  MATRIX *x,  MATRIX *m){
  int i;
  double t1;
  char outfile[2048];
  FILE *fp;

  printf("  Decay data: %s.decay\n", head->outbase);
  if ((strcmp(head->outbase, "-")==0)) fp = stdout;
  else { 
    sprintf(outfile, "%s.decay", head->outbase);
    fp = fopen(outfile, "w");
  }
  OutputHeader(fp, head);
  fprintf(fp, "# Signal normalized coeff.: %g\n", d->mmax);
  fprintf(fp, "# %13s %15s %15s %15s\n",
	  "t1", "Iexp", "Isim", "Iexp - Isim");
  for (i=0; i<d->u.m; i++){
      t1 = d->u.a[i];
      fprintf(fp, "%15.8e %15.8e %15.8e %15.8e\n",
      	      t1, d->m.a[i], m->a[i], d->m.a[i] - m->a[i]);
  }
  fclose(fp);
  fflush(fp);
}

/* ####################################################################### */
/* output for Scurve */
/* ####################################################################### */
void OutputScurve(HEAD *head){
  int i, sw=0;
  double delta;
  char outfile[2048];
  FILE *fp;

  printf("  Scurve data: %s.scurve\n", head->outbase);
  if ((strcmp(head->outbase, "-")==0)) fp = stdout;
  else { 
    sprintf(outfile, "%s.scurve", head->outbase);
    fp = fopen(outfile, "w");
  }
  fprintf(fp, "# %1s %15s %15s %15s %15s %15s\n",
	  "n", "alpha", "sigma", "log10(alpha)", "log10(sigma)", "grad");
  for (i=0; i<head->alpha_num; i++){
    if (i!=0 && (head->alpha_list[i] < head->alpha_list[i-1]) && sw==0){
      fprintf(fp, "\n\n");
      sw = 1;
    }
    fprintf(fp, "%3d %15.8e %15.8e %15.8e %15.8e", 
	    i+1, head->alpha_list[i], head->sigma_list[i],
	    log10(head->alpha_list[i]), log10(head->sigma_list[i]));
    if (i!=(head->alpha_num-1)) {
      delta = (log10(head->sigma_list[i+1]) - log10(head->sigma_list[i]))/
              (log10(head->alpha_list[i+1]) - log10(head->alpha_list[i])); 
      fprintf(fp, " %15.8e\n", delta);
    }
    else  fprintf(fp, "\n");
  }
}

