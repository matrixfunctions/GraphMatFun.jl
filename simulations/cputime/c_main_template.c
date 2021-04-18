#include<stdio.h>
#include<time.h>
#define NOF_SAMPLES 3

int main (int argc, char *argv[]){
  double *A;
  double *E;
  int n,i,j;
  //int n0;
  //clock_t tv1, tv2;
  //double time;

  long seconds,nanoseconds;
  double elapsed;
  struct timespec begin, end;

  //double timingv[NOF_SAMPLES];

  double sum=0;
  clock_gettime(CLOCK_REALTIME, &begin);


  n=2000;
  // n0=10;
  A = malloc(n*n*sizeof(*A));
  E = malloc(n*n*sizeof(*A));
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      A[i+n*j]=0.0;
      if (abs(i-j)<4)

	A[i+n*j]=2.5/8.0;
      if (i==j)
	A[i+n*j]=5.0/8.0;

    }
  }
//  printf("\n");
//  for(i=0; i<n0; i++){
//    for(j=0; j<n0; j++){
//      printf("%f ",A[i+n*j]);
//    }
//    printf("\n");
//  }


  printf("* GRAPH: FUNNAME\n    TIMINGS:  ");
  for (i=0;i<NOF_SAMPLES;i++){
    //tv1=clock();
    clock_gettime(CLOCK_REALTIME, &begin);
    FUNNAME(A, n, E);
    clock_gettime(CLOCK_REALTIME, &end);

    //tv2=clock();
    //time = (tv2 - tv1)/(CLOCKS_PER_SEC);
    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed=seconds + nanoseconds*1e-9;

    printf("%6.3f ",elapsed);
    //timingv[i] = elapsed;
    sum += elapsed;


  }
  printf("\n");
  printf("    MEAN: %.3f\n",sum/ NOF_SAMPLES);
//
//
//  printf("\n");
//  for(i=0; i<n0; i++){
//    for(j=0; j<n0; j++){
//      printf("%f ",E[i+n*j]);
//    }
//    printf("\n");
//  }
  return 0;
//
}
