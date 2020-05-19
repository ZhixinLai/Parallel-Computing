/*
steps of compile and running the code:
gcc gaussian.c -o gaussian -lpthread -lm
./gaussian dim_n p
*/

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <pthread.h>

#define BILLION 1000000000L
#define bool char
#define true 1
#define false 0

struct data_struct {
  int start, end;
  int row;
  int thr_id;
};

double (*A);
double (*B);
double *result;
int maxIndex = 0;
int dim_n = 0;
double maxValue = 0;
pthread_mutex_t mutexmax;
pthread_barrier_t my_barrier;


void *GetMaxValue(void *thr_arg) {

  struct data_struct *mydata;
  mydata = (struct data_struct *) thr_arg;
  int r = mydata->row;
  int s = mydata->start;
  int e = mydata->end;
  int p = 0, index = s;
  double mymax = A[s*dim_n+r];
  for(p = s; p < e; p++) {
    if(mymax < A[p*dim_n+r]) {
      mymax = A[p*dim_n+r];
      index = p;
    }
  }
  pthread_mutex_lock (& mutexmax);
  if(mymax > maxValue) { 
    maxIndex = index;
    maxValue = mymax;
  }
  pthread_mutex_unlock (& mutexmax); 
}


void *GetUpTriangle(void *thr_arg) {

  struct data_struct *mydata;
  mydata = (struct data_struct *) thr_arg;
  int r = mydata->row;
  int s = mydata->start;
  int e = mydata->end;
  int p = 0, q = 0;
  for(p = s; p < e; p++) {
    double rate = A[p*dim_n+r] / A[r*dim_n+r];
    for(q = r; q < dim_n; q++) {
      A[p*dim_n+q] -= rate * A[r*dim_n+q];
    }
    B[p] -= rate * B[r]; 
  } 
}

void *GetDiagonal(void *thr_arg) {

  struct data_struct *mydata;
  mydata = (struct data_struct *) thr_arg;
  int r = mydata->row;
  int s = mydata->start;
  int e = mydata->end;
  int p = 0;
  for(p = s; p < e; p++) {
    double rate = A[p*dim_n+r] / A[r*dim_n+r];
    A[p*dim_n+r] -= rate * A[r*dim_n+r];
    B[p] -= rate * B[r];
  } 
}

void *GetResult(void *thr_arg) {

  struct data_struct *mydata;
  mydata = (struct data_struct *) thr_arg;
  int s = mydata->start;
  int e = mydata->end;
  int p = 0;
  for(p = s; p < e; p++) {
    result[p] = B[p] / A[p*dim_n+p];
  } 
}


int main(int argc, char *argv[]) 
{
	if (argc != 3){
    printf("wrong number of input arguments\n");
    return -1;
	}
  int i, j, k, r, c, thr_p, rc, step, thr_real;
  bool flag = false;
	double drand48(), diff;
	void srand48();
  void *status;
  struct timespec start, end;
  
  srand48(1);
	dim_n = atoi(argv[1]);
  thr_p = atoi(argv[2]);
  pthread_t thread[thr_p];
  pthread_mutex_init(&mutexmax, NULL);
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_barrier_init(&my_barrier,NULL,2);
  
  // set up the data in the heap
	double* A_temp = (double*)malloc(dim_n * dim_n * sizeof(double));
  double* B_temp = (double*)malloc(dim_n * sizeof(double));
  double* result_temp = (double*)malloc(dim_n * sizeof(double));
  double tmp;
  struct data_struct datas[thr_p];
	// populate matrices
	for (i = 0; i < dim_n; i++){
		for (j = 0; j < dim_n; j++){
			A_temp[i*dim_n+j] = drand48();
		}
    B_temp[i] = drand48();
    result_temp[i] = drand48();
	}
  A = A_temp;
  B = B_temp;
  result = result_temp;

	clock_gettime(CLOCK_MONOTONIC, &start);
 
  /* step1 triangularization */
  for(i = 0; i < dim_n; i++){
  
    // get the max value and index
    if((dim_n - i) % thr_p == 0) step = (dim_n - i) / thr_p;
    else step = (dim_n - i) / thr_p + 1;
    maxValue = A[i*dim_n+i];
    maxIndex = i;
    thr_real = thr_p;
    flag = false;
    for(k = 0; k < thr_p; k++){
      datas[k].row = i;
      datas[k].start = step * k + i;
      datas[k].end = step * (k + 1) + i;
      datas[k].thr_id = k;
      if(datas[k].end >= dim_n) {
        datas[k].end = dim_n;
        thr_real = k + 1;
        flag = true;
      }
      rc = pthread_create(&thread[k], &attr, GetMaxValue, (void *) &datas[k]);
      if(flag) break;
    }
    for(k = 0; k < thr_real; k++) {
      rc = pthread_join(thread[k], &status);
    }
    
    // swap the line i with the maximum one
    for(k = i; k < dim_n; k++) {
      tmp = A[maxIndex*dim_n+k];
      A[maxIndex*dim_n+k] = A[i*dim_n+k];
      A[i*dim_n+k] = tmp;
    }
    tmp = B[i];
    B[i] = B[maxIndex];
    B[maxIndex] = tmp;

    // update a(i+1:n,i+1:n) and b(i:n)
    if((dim_n - i - 1) % thr_p == 0) step = (dim_n - i - 1) / thr_p;
    else step = (dim_n - i - 1) / thr_p + 1;
    thr_real = thr_p;
    flag = false;
    for(k = 0; k < thr_p; k++){
      datas[k].row = i;
      datas[k].start = step * k + i + 1;
      datas[k].end = step * (k + 1) + i + 1;
      datas[k].thr_id = k;
      
      if(datas[k].end >= dim_n) {
        datas[k].end = dim_n;
        thr_real = k + 1;
        flag = true;
      }
      rc = pthread_create(&thread[k], &attr, GetUpTriangle, (void *) &datas[k]);
      if(flag) break;
    }
    for(k = 0; k < thr_real; k++) {
      rc = pthread_join(thread[k], &status);
    }  
  }


  /* step2 Diagonalization */
  for(i = dim_n - 1; i >= 0; i--){
    
    if(i % thr_p == 0) step = (i) / thr_p;
    else step = (i) / thr_p + 1;
    thr_real = thr_p;
    flag = false;
    for(k = 0; k < thr_p; k++){
      datas[k].row = i;
      datas[k].start = step * k;
      datas[k].end = step * (k + 1);
      datas[k].thr_id = k;

      if(datas[k].end >= i) {
        datas[k].end = i;
        thr_real = k + 1;
        flag = true;
      }
      rc = pthread_create(&thread[k], &attr, GetDiagonal, (void *) &datas[k]);
      if(flag) break;
    }
    for(k = 0; k < thr_real; k++) {
      rc = pthread_join(thread[k], &status);
    }  
  }

  /* step3 Get the result */
  if(dim_n % thr_p == 0) step = (dim_n) / thr_p;
  else step = (dim_n) / thr_p + 1;
  thr_real = thr_p;
  flag = false;
  for(k = 0; k < thr_p; k++){
    datas[k].start = step * k;
    datas[k].end = step * (k + 1);
    datas[k].thr_id = k;

    if(datas[k].end >= dim_n) {
      datas[k].end = dim_n;
      thr_real = k + 1;
      flag = true;
    }
    rc = pthread_create(&thread[k], &attr, GetResult, (void *) &datas[k]);
    if(flag) break;
  }
  for(k = 0; k < thr_real; k++) {
    rc = pthread_join(thread[k], &status);
  }

  
  pthread_attr_destroy(&attr);
  pthread_barrier_destroy(&my_barrier);
	clock_gettime(CLOCK_MONOTONIC,&end);
	diff = BILLION * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
  printf("the dim_n is %d, the number of thread is %d, the execution time is %lf nanoseconds\n", dim_n, thr_p, diff);
  double different = 0;
  
  /* step4 check the answer */
  for(r = 0; r < dim_n; r++) {
    tmp = 0; 
    for(c = 0; c < dim_n; c++) {
        tmp += A[r*dim_n+c] * result[c];
    }
    different += pow((B[r] - tmp), 2); 
  }
  different = pow(different, 0.5);
  printf("the difference between the real result and prediction is %f\n", different);
  
	return 0;
}
