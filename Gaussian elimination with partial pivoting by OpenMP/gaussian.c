/*
steps of compile and unning the code:
gcc gaussian.c -o gaussian -fopenmp -lm
gaussian dim_n p
*/

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <omp.h>


#define BILLION 1000000000L
#define bool char

double rate;
int maxIndex = 0;
int myindex = 0;
int dim_n = 0;
double maxValue = 0;
double mymax = 0;
double* A ;
double* A_test;
double* B;
double* B_test;
double* result;
double* RR;
double tmp;
int i, j, k, r, c, p, q, thr_p, rc, step;

void getMaxwithOpenMP(int i){
  #pragma omp parallel shared(maxValue, maxIndex) private(mymax, myindex) 
  {
    mymax = A[i*dim_n+i];
    myindex = i;
    #pragma omp for schedule(dynamic)
      for(k = i+1; k < dim_n; k++) {
        if(mymax < A[k*dim_n+i]) {
          mymax = A[k*dim_n+i];
          myindex = k;
        }
      }
    #pragma omp critical(name) 
    {
      if(mymax > maxValue) { 
        maxIndex = myindex;
        maxValue = mymax;
      }
    }
  }
}

void getMaxwithoutOpenMP(int i){
  mymax = A[i*dim_n+i];
  myindex = i;
  for(k = i+1; k < dim_n; k++) {
    if(mymax < A[k*dim_n+i]) {
      mymax = A[k*dim_n+i];
      myindex = k;
    }
  }
  if(mymax > maxValue) { 
    maxIndex = myindex;
    maxValue = mymax;
  }
}

void swapwithMax(int i){
  for(k = i; k < dim_n; k++) {
    tmp = A[maxIndex*dim_n+k];
    A[maxIndex*dim_n+k] = A[i*dim_n+k];
    A[i*dim_n+k] = tmp;
  }
  for(k = 0; k < dim_n; k++) {
    tmp = B[maxIndex*dim_n+k];
    B[maxIndex*dim_n+k] = B[i*dim_n+k];
    B[i*dim_n+k] = tmp;
  }
}

void getUpTrianglewithOpenMP(int i){
  #pragma omp parallel for private(rate, p, q) schedule(dynamic)
    for(p = i+1; p < dim_n; p++) {
      rate = A[p*dim_n+i] / A[i*dim_n+i];
      for(q = i; q < dim_n; q++) {
        A[p*dim_n+q] -= rate * A[i*dim_n+q];
      }
      for(q = 0; q < dim_n; q++) {
        B[p*dim_n+q] -= rate * B[i*dim_n+q];
      }
    } 
}

void getUpTrianglewithoutOpenMP(int i){
  for(p = i+1; p < dim_n; p++) {
    rate = A[p*dim_n+i] / A[i*dim_n+i];
    for(q = i; q < dim_n; q++) {
      A[p*dim_n+q] -= rate * A[i*dim_n+q];
    }
    for(q = 0; q < dim_n; q++) {
      B[p*dim_n+q] -= rate * B[i*dim_n+q];
    }
  } 
}
    
void backSubbstitutionwithOpenMP(){
  for(i = dim_n - 1; i >= 0; i--){
    #pragma omp parallel for private(p) schedule(dynamic)
    for(p = 0; p < dim_n; p++){
      result[i*dim_n+p] = B[i*dim_n+p] / A[i*dim_n+i];
    }
    #pragma omp parallel for private(j, p) schedule(dynamic)
    for(j = i-1; j >= 0; j--) {
      for(p = 0; p < dim_n; p++){
        B[j*dim_n+p] -= A[j*dim_n+i] * result[i*dim_n+p];
      }
    }
  }
}

void backSubbstitutionwithoutOpenMP(){
  for(i = dim_n - 1; i >= 0; i--){
    for(p = 0; p < dim_n; p++){
      result[i*dim_n+p] = B[i*dim_n+p] / A[i*dim_n+i];
    }
    for(j = i-1; j >= 0; j--) {
      for(p = 0; p < dim_n; p++){
        B[j*dim_n+p] -= A[j*dim_n+i] * result[i*dim_n+p];
      }
    }
  }
}
  
double checkAnswer(){
  double different = 0, nR = 0, nA = 0, nX = 0;
  for(r = 0; r < dim_n; r++) {
    for(c = 0; c < dim_n; c++) {
      for(p = 0; p < dim_n; p++) {
        RR[r*dim_n+c] += A_test[r*dim_n+p] * result[p*dim_n+c];
      }
      RR[r*dim_n+c] -= B_test[r*dim_n+c];
      nR += pow(RR[r*dim_n+c], 2);
      nX += pow(result[r*dim_n+c], 2);
      nA += pow(A_test[r*dim_n+c], 2); 
    }
  }
  different = pow(nR, 0.5) / (pow(nX, 0.5) * pow(nA, 0.5));
  return different;
}

int main(int argc, char *argv[]) 
{
	if (argc != 3){
    printf("wrong number of input arguments\n");
    return -1;
	}
 
	double drand48(), diff;
	void srand48();
  void *status;
  struct timespec start, end;
  
  srand48(1);
	dim_n = atoi(argv[1]);
  thr_p = atoi(argv[2]);
  omp_set_num_threads(thr_p);
  
  // set up the data in the heap
	A = (double*)malloc(dim_n * dim_n * sizeof(double));
  A_test = (double*)malloc(dim_n * dim_n * sizeof(double));
  B = (double*)malloc(dim_n * dim_n * sizeof(double));
  B_test = (double*)malloc(dim_n * dim_n * sizeof(double));
  result = (double*)malloc(dim_n * dim_n * sizeof(double));
  RR = (double*)malloc(dim_n * dim_n * sizeof(double));
  
  
	// populate matrices
	for (i = 0; i < dim_n; i++){
		for (j = 0; j < dim_n; j++){
			A[i*dim_n+j] = drand48();
      A_test[i*dim_n+j] = A[i*dim_n+j];
      if(i == j) B[i*dim_n+j] = 1;
      else B[i*dim_n+j] = 0;
      B_test[i*dim_n+j] = B[i*dim_n+j];
      RR[i*dim_n+j] = 0;
      result[i*dim_n+j] = drand48();
		}
	}

	clock_gettime(CLOCK_MONOTONIC, &start);
 
  /* step1 triangularization */
  for(i = 0; i < dim_n; i++){
  
    // get the max value and index
    maxValue = A[i*dim_n+i];
    maxIndex = i;
    
    if(thr_p == 0) getMaxwithoutOpenMP(i);
    else getMaxwithOpenMP(i);
    
    // swap the line i with the maximum one
    swapwithMax(i);
        
    // update a(i+1:n,i+1:n) and b(i+1:n, :)
    if(thr_p == 0) getUpTrianglewithoutOpenMP(i);
    else getUpTrianglewithOpenMP(i);    
  }
  
  /* step2 backsubbstitution */
  if(thr_p == 0) backSubbstitutionwithoutOpenMP();
  else backSubbstitutionwithOpenMP();

	clock_gettime(CLOCK_MONOTONIC,&end);
	diff = BILLION * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
  printf("the dim_n is %d, the number of thread is %d, the execution time is %e nanoseconds\n", dim_n, thr_p, diff);
  // printf("%d,%d,%e\n", dim_n, thr_p, diff);
  
  /* step3 check the answer */
  double different = checkAnswer();
  if(different < pow(10, -12)) printf("the difference between the real result and prediction is less than 10^-12\n");

  return 0;

}
