#include <unistd.h>
#include <getopt.h>
#include "ilt.h"

void Error(){
  printf("Usage: osilap [options] infile.int\n\n");
  printf("Options:\n");
  printf(" -i                      confirm calculation condition\n");
  printf(" -h --help               show this help message\n");
  printf("    --dimension=%%d       one- or two-dimension [2]\n");
  printf("    --alpha=%%f           initial value of alpha (only BRD) [0.1]\n");
  printf("    --sigma=%%f           standard error, 0=auto [0.0]\n");
  printf("    --sfactor=%%f         sfactor for SVD [1e-4]\n");
  printf("    --output=%%s          basename of output file [basename]\n");
  printf("    --TXspace=%%s         TX space (log | linear) [log]\n");
  printf("    --TXtype=%%s          TX data type (T1 | T1inf | T1p | T2 | D) [T1]\n");
  printf("    --TXmin=%%f           TX minimum [0.0001]\n");
  printf("    --TXmax=%%f           TX maximum [10]\n");
  printf("    --TXnum=%%d           TX points [50]\n");
  printf("    --TYspace=%%s         TY space (log | linear) [log]\n");
  printf("    --TYtype=%%s          TY data type (T1 | T2 | D) [T2]\n");
  printf("    --TYmin=%%f           TY minimum [0.0001]\n");
  printf("    --TYmax=%%f           TY maximum [10]\n");
  printf("    --TYnum=%%d           TY points [50]\n");
  printf("    --alpha_loop_max=%%d  maximum iteration for alpha [5000]\n");
  printf("    --alpha_tol=%%f       alpha tolerance [1e-3]\n");
  printf("    --newton_loop_max=%%d maximum iteration for newton method [5000]\n");
  printf("    --newton_tol=%%f      newton tolerance [1e-8]\n");
  printf("    --Scurve=%%f          target for dlog10(chi)/dlog10(alpha) [0.0]\n");
  printf("    --Smin=%%f            search grid for s-curve [1e+6]\n");
  printf("    --Smax=%%f            search grid for s-curve [1e-6]\n");
  printf("    --Snum=%%d            search grid number for s-curve [20]\n");
  printf("    --Atype=%%d           calculation method of alpha for BRD routine [0]\n");
  exit(1);
}

void Confirm(HEAD *head, INDATA *d){
  int ch;
  printf("\n");
  printf("-------- Solving condition -------\n");
  printf(" dimension:       %d\n", head->dimension);
  printf(" sfactor:         %g\n", head->sfactor);
  printf(" initial alpha:   %g\n", head->alpha0);
  printf(" sigma:           %g\n", head->sigma0);
  printf(" alpha max Loop:  %d\n", head->alpha_loop_max);
  printf(" alpha tol:       %g\n", head->alpha_tol);
  printf(" Newton max Loop: %d\n", head->newton_loop_max);
  printf(" Newton tol:      %g\n", head->newton_tol);
  printf(" A_opt method:    %d\n", head->Atype);
  printf("--------- TX information ---------\n");
  printf(" Grid space:      %s\n", head->TXspace);
  printf(" Data type:       %s\n", head->TXtype);
  printf(" (TXmin, TXmax, TXpoints) =");
  printf(" (%g, %g, %d)\n", head->TXmin, head->TXmax, head->TXnum);
  if (head->dimension>1){
    printf("--------- TY information ---------\n");
    printf(" Grid space:      %s\n", head->TYspace);
    printf(" Data type:       %s\n", head->TYtype);
    printf(" (TYmin, TYmax, TYpoints) =");
    printf(" (%g, %g, %d)\n", head->TYmin, head->TYmax, head->TYnum);
  }
  if (head->s_taget != 0){
  printf("------ S-curve configuration -----\n");
  printf(" Alpha grid space: \n");
  printf(" (start, end, points) =");
  printf(" (%g, %g, %d)\n", head->alpha_start, head->alpha_end, head->alpha_num);
  }
  printf("----------- Input data -----------\n");
  printf(" Input file:      %s\n", head->infile);
  printf(" Input data size:");
  if (head->dimension == 1){
    printf(" %d\n", d->u.m);
  }
  if (head->dimension == 2){
    printf(" (X, Y)=(%d, %d)\n", d->u.m, d->v.m);
  }
  printf(" Output basename: %s\n", head->outbase);
  printf("\n");

  fflush(stdout);
  if (head->confirm == 1){
    printf("Start calculation? [y/n]: ");
    while ((ch = getchar()) != EOF) {        /* 文字の入力 */
      if (ch == 'y') break;   
      if (ch == 'n') exit(1);
    }
  }
  fflush(stdin);
}

/* ####################################################################### */
/* ベース名を決める */
/* ####################################################################### */
void GetBaseName(HEAD *head)
{
  int nLen;
  char *res = NULL;
  char szPath[1024];

  strcpy(szPath, head->infile);
  nLen = sizeof (szPath);
  memset(head->outbase, 0, sizeof(head->outbase));
  if ((nLen > 0) && (nLen <= 1024)) {
    res = strrchr( szPath, '.' );
    if(res != NULL)
      *res = '\0';
  }
  strcpy(head->outbase, szPath);
}


int SetArgment(int argc, char * argv[], HEAD *head){
  int opt, option_index;
  struct option long_options[] = {
    {"help",  no_argument,  NULL, 'h'},
    {"confirm",  no_argument,  NULL, 'i'},
    {"DX",  required_argument,  NULL, 'A'},
    {"DY",  required_argument,  NULL, 'B'},
    {"TXtype",  required_argument,  NULL, 'C'},
    {"TXspace",  required_argument,  NULL, 'D'},
    {"TXmin",  required_argument,  NULL, 'E'},
    {"TXmax",  required_argument,  NULL, 'F'},
    {"TXnum",  required_argument,  NULL, 'G'},
    {"TYtype",  required_argument,  NULL, 'H'},
    {"TYspace",  required_argument,  NULL, 'I'},
    {"TYmin",  required_argument,  NULL, 'J'},
    {"TYmax",  required_argument,  NULL, 'K'},
    {"TYnum",  required_argument,  NULL, 'L'},
    {"output",  required_argument,  NULL, 'M'},
    {"sfactor",  required_argument,  NULL, 'N'},
    {"dimension",  required_argument,  NULL, 'O'},
    {"alpha",  required_argument,  NULL, 'P'},
    {"alpha_loop_max",  required_argument,  NULL, 'Q'},
    {"alpha_tol",  required_argument,  NULL, 'R'},
    {"newton_loop_max",  required_argument,  NULL, 'S'},
    {"newton_tol",  required_argument,  NULL, 'T'},
    {"sigma",  required_argument,  NULL, 'U'},
    {"Scurve",  required_argument,  NULL, 'V'},
    {"Smin",  required_argument,  NULL, 'W'},
    {"Smax",  required_argument,  NULL, 'X'},
    {"Snum",  required_argument,  NULL, 'Y'},
    {"Atype",  required_argument,  NULL, 'Z'},
    {0, 0, 0, 0}// 配列の最後はすべて0で埋める
  };
    
  //opterr = 0;/* エラーメッセージを非表示にする */
  if (argc == 1) Error();
    
  //while((opt = getopt(argc, argv, "mes:")) != -1){
  while((opt = getopt_long(argc, argv, "ih", 
			   long_options, &option_index)) != -1){
    switch(opt){
    case 'i':
      head->confirm = 1;
      break;
    case 'A': head->DXnum = atoi(optarg);  break;
    case 'B': head->DYnum = atoi(optarg);  break;
    case 'C': strcpy(head->TXtype, optarg);  break;
    case 'D': strcpy(head->TXspace, optarg); break;
    case 'E': head->TXmin = atof(optarg); break;
    case 'F': head->TXmax = atof(optarg); break;
    case 'G': head->TXnum = atoi(optarg); break;
    case 'H': strcpy(head->TYtype, optarg);  break;
    case 'I': strcpy(head->TYspace, optarg); break;
    case 'J': head->TYmin = atof(optarg); break;
    case 'K': head->TYmax = atof(optarg); break;
    case 'L': head->TYnum = atoi(optarg); break;
    case 'M': strcpy(head->outbase, optarg); break;
    case 'N': head->sfactor = atof(optarg);  break;
    case 'O': head->dimension = atoi(optarg);  break;
    case 'P': head->alpha0 =  atof(optarg);  break;
    case 'Q': head->alpha_loop_max = atoi(optarg);  break;
    case 'R': head->alpha_tol = atoi(optarg); break;
    case 'S': head->newton_loop_max = atoi(optarg); break;
    case 'T': head->newton_tol = atof(optarg);  break;
    case 'U': head->sigma0 = atof(optarg);   break;
    case 'V': head->s_taget = atof(optarg);  break;
    case 'W': head->alpha_start = atof(optarg);  break;
    case 'X': head->alpha_end = atof(optarg);  break;
    case 'Y': head->alpha_num = atoi(optarg);  break;
    case 'Z': head->Atype = atoi(optarg);  break;
      // 解析できないオプションが見つかった場合は「?」を返す
      // オプション引数が不足している場合も「?」を返す
    case '?':
      Error();
    }
  }
  strcpy(head->infile, argv[optind]);
  if (strcmp(head->outbase, "")==0) GetBaseName(head);
  /* if (head->dimension == 1 && head->TXnum < 100){ */
  /*   head->TXnum = 200; */
  /*   strcpy(head->TXtype, "T2"); */
  /* } */
  return 0; 
}
