#include "ilt.h"

void SetDataSize(HEAD *head, INDATA *d){
  char buf[1024];
  char *tp;

  if((head->infp = fopen(head->infile, "r"))==NULL){
    printf("%s is not found.\n", head->infile);
    exit(1);
  }

  tp = fgets(buf, sizeof(buf), head->infp);
  if (tp == NULL){
    printf("Read error. %s\n", head->infile);
    exit(1);
  }
  rewind(head->infp);

  if (head->DXnum == 0 && head->DYnum == 0){
    /* スペース.を区切りに文字列を抽出    */
    tp = strtok(buf, "X:/Y:" );
    tp = strtok(NULL, "X:/Y:" );
    if (tp == NULL){
      printf("X Dimesion for input data error!\n");
      exit(1);
    }
    head->DXnum = atoi(tp);  /* tau1軸がu */;
    tp = strtok(NULL, "X:/Y:" );
    if (tp == NULL){
      printf("Y Dimesion for input data error!\n");
      exit(1);
    }
    head->DYnum = atoi(tp);  /* tau2軸がv */
  }
  if (head->DXnum == 0 || head->DYnum == 0){
    printf("XY Dimesion for input data error!\n");
    exit(1);
  }
  d->u.m =  head->DXnum;
  d->v.m =  head->DYnum;
  d->u.n =  1;
  d->v.n =  1;
  d->u.a = (double*)malloc(sizeof(double) * d->u.m );
  d->v.a = (double*)malloc(sizeof(double) * d->v.m);
  d->m.a = (double*)malloc(sizeof(double) * d->u.m * d->v.m);
  d->m.m = d->u.m;  /* T1軸が行 */
  d->m.n = d->v.m;  /* T2軸が列 */
}


/* ####################################################################### */
/* 1Dの測定データを読む */
/* ####################################################################### */
void ReadData1D(HEAD *head, INDATA *d){
  char buf[1024];
  int pnt = 0, col;

  if ((strcmp(head->infile, "-")==0)) head->infp = stdin;
  else if((head->infp = fopen(head->infile, "r"))==NULL){
    printf("%s is not found.\n", head->infile);
    exit(1);
  }
  if (head->DXnum == 0){
    while (fgets(buf, sizeof(buf), head->infp) != NULL){
      if(strstr(buf, "#") == NULL && buf[0] != '\n'){
	pnt++;
      }
    }
    head->DXnum = pnt;
    rewind(head->infp);
  }
  d->u.m =  head->DXnum;
  d->u.n =  1;
  d->v.m =  0;
  d->v.n =  0;
  d->u.a = (double*)malloc(sizeof(double) * head->DXnum);
  d->m.a = (double*)malloc(sizeof(double) * head->DXnum);
  d->m.m = d->u.m; 
  d->m.n = 1;  
  head->s1 = head->DXnum;
  head->s2 = 0;

  pnt = 0;
  while (fgets(buf, sizeof(buf), head->infp) != NULL){
    if(strstr(buf, "#") == NULL && buf[0] != '\n' && buf[0] != '\r'){
      col = sscanf(buf, "%lf %lf", &d->u.a[pnt], &d->m.a[pnt]);
      if (col != 2){
	printf("read error. line numeber=%d (%s)\n", pnt, head->infile);
	exit(1);
      }
      /* d->u.a[pnt] = d->u.a[pnt] * 1e-3;   /\* s->ms *\/ */

      pnt++;
      if (pnt == head->DXnum) break;
    }
  }

  fclose(head->infp);

  /* 最大値で規格化 */
  /* d->mmax = fabs(d->m.a[idamax_(&m, d->m.a, &incr)-1]); */
  /* alpha = 1.0/d->mmax; */
  /* dscal_(&m, &alpha, d->m.a, &incr); */
}

/* ####################################################################### */
/* 測定データの読み込み  */
/* ####################################################################### */
void ReadData(HEAD *head, INDATA *d){
  char buf[1024];
  int i=0, j=0, p, col;
  double ta, tb, mi;
  div_t temp;
  
  while (fgets(buf, sizeof(buf), head->infp) != NULL){
    if(strstr(buf, "#") == NULL && buf[0] != '\n' && buf[0] != '\r'){
      col = sscanf(buf, "%lf %lf %lf", &ta, &tb, &mi);
      if (col != 3){
	printf("read error. line:%d (%s)\n", i, head->infile);
	exit(1);
      }
      /* X軸一定でY軸変化のためM行列を転置 */
      temp = div(i, d->m.n);
      p = temp.rem * d->m.m + temp.quot;
      /* Mの並びがT1変化が先なら p=iにする */
      d->m.a[p] = mi;
      if (i < d->v.m) d->v.a[i] = tb;
      if ((i % d->v.m)==0) {
	d->u.a[j] = ta;
	j++;
      }
      i++;
    }
  }

  /* mn = d->m.m * d->m.n; */
  /* d->mmax = fabs(d->m.a[idamax_(&mn, d->m.a, &incr)-1]); */
  /* alpha = 1.0/d->mmax; */
  /* dscal_(&mn, &alpha, d->m.a, &incr); */
}
