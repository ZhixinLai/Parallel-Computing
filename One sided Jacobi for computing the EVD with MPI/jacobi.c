/*
steps of compile and unning the code:
(1) no MPI
mpicc jacobi.c -o jacobi -lm
./jacobi 128 10000 -16 0  

*explanation*:     
// 128 is dim_n, 10000 is the max iter number;                                    
// -16 means the tol = 10^-16; 0 means not use MPI.

(2) use MPI

mpicc jacobi.c -o jacobi -lm
sbatch jacobi.sub  

*attention*: 
// in .sub file: mem-per-cpu should be changed according to the tasks-per-node and dim_n
// if dim_n is high, the mem-per-cpu should be high, otherwise memory limit will occur
// mem-per-cpu * tasks-per-node must be less 90000, "otherwise Memory specification can not be satisfied" error will happen
// as for most cases, mem-per-cpu = 4000 / 8000 will be ok; as for dim_n = 1024, mem-per-cpu should be higher
 
*explanation*:    
// in .sh file: ./jacobi 128 10000 -16 1      
// 128 is dim_n, 10000 is the max iter number; 
// -16 means the tol = 10^-16; 1 means use MPI
*/

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>


#define a(x, y) a[x * 2 + y]
#define A(x, y) A[x * dim_n + y]
#define V(x, y) V[x * dim_n + y]
#define ALocal(x, y) ALocal[x * dim_n + y]
#define VLocal(x, y) VLocal[x * dim_n + y]
#define ALocal2(x, y) ALocal2[x * dim_n + y]
#define VLocal2(x, y) VLocal2[x * dim_n + y]
#define ALocalTmp(x, y) ALocalTmp[x * dim_n + y]
#define VLocalTmp(x, y) VLocalTmp[x * dim_n + y]
#define ALocal2Tmp(x, y) ALocal2Tmp[x * dim_n + y]
#define VLocal2Tmp(x, y) VLocal2Tmp[x * dim_n + y]
#define GT(x, y) GT[x * 2 + y]
#define B(x, y) B[x * dim_n + y]
#define Sigma(x, y) Sigma[x * dim_n + y]

double* A, *V, *ALocal, *VLocal, *l_buf_a, *r_buf_a, *l_buf_v, *r_buf_v, *sumLocal;
int P_ROWS, dim_n, max_sweeps, N, num_el, iter_n, sum_gap;
int rank, npes, right, left, row_size, i, j, MPI_flag;
double tol = 0, sum = 10, threshold = 0, threshold_tmp = 0;
struct timespec start, end;


MPI_Datatype row_type;
MPI_Status stats[2];
MPI_Request reqs[2];

void generateALocalVLocal() {

  ALocal = (double *) malloc(num_el*sizeof(double));
  VLocal = (double *) malloc(num_el*sizeof(double));
  l_buf_a = (double *) calloc(dim_n, sizeof(double));
  r_buf_a = (double *) calloc(dim_n, sizeof(double));
  l_buf_v = (double *) calloc(dim_n, sizeof(double));
  r_buf_v = (double *) calloc(dim_n, sizeof(double));
  if(rank == 0 || MPI_flag == 0) {
    A = (double *) malloc(N*sizeof(double));
    V = (double *) malloc(N*sizeof(double));

    for(i = 0; i<dim_n; i++){ 
      for(j=0;j<dim_n;j++) {
        A[i*dim_n+j] = drand48();
        if(i == j) {
          V[i*dim_n+j] = 1;
        }
        else V[i*dim_n+j] = 0;
        // printf("%f ",A[i*dim_n+j]);
      }
      // printf("\n");
    }  
  }

  // Scatter to each of npes processes if use MPI 
  if(MPI_flag == 1) {
    MPI_Scatter(A, num_el, MPI_DOUBLE, ALocal, num_el, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    MPI_Scatter(V, num_el, MPI_DOUBLE, VLocal, num_el, MPI_DOUBLE, 0, MPI_COMM_WORLD);   
  }
  // ALocal = A if not use MPI 
  else {
    for(i = 0; i<dim_n; i++){ 
      for(j=0;j<dim_n;j++) {
          ALocal(i,j) = A(i,j);
          VLocal(i,j) = V(i,j);
      }
    }  
  }  
}

void rotationComputation(double* ALocal2, double* VLocal2) {
  
  double* a = (double*)malloc(2 * 2 * sizeof(double));
  int i = 0, j = 0, k = 0;
  
  for(i = 0; i < 2; i ++) {
    for(j = 0; j < 2; j++) {
      a(i, j) = 0;
      for(k = 0; k < dim_n; k++) {
        a(i, j) += ALocal2(i, k) * ALocal2(j, k);
      }
    }
  } 
  
  double tau = (a(1,1)-a(0,0))/(2*a(0,1));
  double flag = 0;
  if(tau > 0) flag = 1;
  else flag = -1; 
  double t = 1/(tau+flag*pow(1+tau*tau,0.5));
  double c = 1/pow(1+t*t,0.5), s = c*t;
  double* GT = (double*)malloc(2 * 2 * sizeof(double));
  GT(0, 0) = c; GT(0, 1) = -s; GT(1, 0) = s; GT(1, 1) = c;
  double* ALocal2Tmp = (double*)malloc(2 * dim_n * sizeof(double));
  double* VLocal2Tmp = (double*)malloc(2 * dim_n * sizeof(double));
  
  for(i = 0; i < 2 ; i++) {
    for(j = 0; j < dim_n; j++) {
      ALocal2Tmp(i, j) = 0;
      VLocal2Tmp(i, j) = 0;
      for(k = 0; k < 2; k++) {
        ALocal2Tmp(i, j) +=  GT(i, k) * ALocal2(k, j);
        VLocal2Tmp(i, j) +=  GT(i, k) * VLocal2(k, j);
      }
    }
  }

  for(i = 0; i < 2 ; i++) {
    for(j = 0; j < dim_n; j++) {
      ALocal2(i, j) = ALocal2Tmp(i, j);
      VLocal2(i, j) = VLocal2Tmp(i, j);
    }
  }
  sumLocal[0] += 2*a(0,1)*a(0,1);
  sumLocal[1] += a(0,0)*a(0,0)+a(1,1)*a(1,1);
}


void switchRows(){
  // [1] use MPI and npes > 1
  if(MPI_flag == 1 && npes != 1){
  // (1) send right, receive from left
  // general case(not PE0 nor npes-1)

    if ((rank>0)&&(rank<npes-1)){
      right = rank+1; left = rank-1;
      MPI_Send(&ALocal[(P_ROWS-2)*dim_n],1,row_type,right,0,MPI_COMM_WORLD);
      MPI_Recv(r_buf_a,1,row_type,left,0,MPI_COMM_WORLD,&stats[0]);
      MPI_Send(&VLocal[(P_ROWS-2)*dim_n],1,row_type,right,1,MPI_COMM_WORLD);
      MPI_Recv(r_buf_v,1,row_type,left,1,MPI_COMM_WORLD,&stats[1]);
    }
  // boundary case for PE0, leftmost PE

    if (rank==0){
      right = rank+1; left = rank;
      if(P_ROWS == 2) {
        MPI_Send(&ALocal[(P_ROWS-1)*dim_n],1,row_type,right,0, MPI_COMM_WORLD);
        MPI_Send(&VLocal[(P_ROWS-1)*dim_n],1,row_type,right,1, MPI_COMM_WORLD);
      }
      else {
        MPI_Send(&ALocal[(P_ROWS-2)*dim_n],1,row_type,right,0, MPI_COMM_WORLD);
        MPI_Send(&VLocal[(P_ROWS-2)*dim_n],1,row_type,right,1, MPI_COMM_WORLD);
      }
    }
  // boundary case for PE(npes-1), rightmost PE

    if (rank == npes-1){
      left = rank-1; right = rank;
      // MPI_Send(&ALocal[(P_ROWS-2)*dim_n],1,row_type,right,0,MPI_COMM_WORLD);
      MPI_Recv(r_buf_a,1,row_type,left,0,MPI_COMM_WORLD,&stats[0]);
      // MPI_Send(&VLocal[(P_ROWS-2)*dim_n],1,row_type,right,1,MPI_COMM_WORLD);
      MPI_Recv(r_buf_v,1,row_type,left,1,MPI_COMM_WORLD,&stats[1]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  
  // (2) send left, receive from right
  // general case(not PE0 nor npes-1)

    if ((rank>0)&&(rank<npes-1)){
      right = rank+1; left = rank-1;
      MPI_Send(&ALocal[dim_n],1,row_type,left,0,MPI_COMM_WORLD);
      MPI_Recv(l_buf_a,1,row_type,right,0,MPI_COMM_WORLD,&stats[0]);
      MPI_Send(&VLocal[dim_n],1,row_type,left,1,MPI_COMM_WORLD);
      MPI_Recv(l_buf_v,1,row_type,right,1,MPI_COMM_WORLD,&stats[1]);
    }
  // boundary case for PE0, leftmost PE
    if (rank==0){
      right = rank+1; left = rank;
      MPI_Recv(l_buf_a,1,row_type,right,0, MPI_COMM_WORLD,&stats[0]);
      MPI_Recv(l_buf_v,1,row_type,right,1, MPI_COMM_WORLD,&stats[1]);
    }
  // boundary case for PE(mpes-1), rightmost PE
    if (rank == npes-1){
      left = rank-1; right = rank;
      MPI_Send(&ALocal[dim_n],1,row_type,left,0,MPI_COMM_WORLD);
      MPI_Send(&VLocal[dim_n],1,row_type,left,1,MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    // store templely the mostright row 8
    double* ALocalRM = (double*)malloc(dim_n * sizeof(double));
    double* VLocalRM = (double*)malloc(dim_n * sizeof(double));
    for(i = 0; i < dim_n; i++){
      ALocalRM[i] = ALocal((P_ROWS-1), i);
      VLocalRM[i] = VLocal((P_ROWS-1), i);
    }

    // exchange within a PE
  
    // as for rank = npes-1, assign the second mostright to mostright 7->8
    if(rank == npes-1) {
      for(j = 0; j < dim_n; j++) {
        ALocal((P_ROWS-1), j) = ALocal((P_ROWS-2), j);
        VLocal((P_ROWS-1), j) = VLocal((P_ROWS-2), j);      
      }
    }
    
    // assign left to right 5->7
    for(i = P_ROWS-2; i > 0; i=i-2) {
      for(j = 0; j < dim_n; j++) {
        if(rank == 0 && i == 2) {
          ALocal(i, j) = ALocal((i-1), j);
          VLocal(i, j) = VLocal((i-1), j);
        }
        else {
          ALocal(i, j) = ALocal((i-2), j);
          VLocal(i, j) = VLocal((i-2), j);
        }
      }
    }
    // assign right to left 6->4
    for(i = 1; i < P_ROWS-2; i=i+2) {
      for(j = 0; j < dim_n; j++) {
        // mostright element
        if(i == P_ROWS-3) {
          ALocal(i, j) = ALocalRM[j];
          VLocal(i, j) = VLocalRM[j];
        }
        else {
          ALocal(i, j) = ALocal((i+2), j);
          VLocal(i, j) = VLocal((i+2), j);
        }
      }
    }
  
  
  // exchange among different PEs
    // from right
    if(rank != 0)
      for(j = 0; j < dim_n; j++) {
        ALocal(0, j) = r_buf_a[j];
        VLocal(0, j) = r_buf_v[j];
      }
    // from left
    if(rank != npes-1)
      for(j = 0; j < dim_n; j++) {
        ALocal((P_ROWS-1), j) = l_buf_a[j];
        VLocal((P_ROWS-1), j) = l_buf_v[j];
      }
  }  
  
  // [2] not use MPI or npes = 1
  if(npes == 1 || MPI_flag == 0) {
  
  // exchange within a PE
    // store the rightmost-1 row temp
    double *ALocalR = (double *) calloc(dim_n, sizeof(double));
    double *VLocalR = (double *) calloc(dim_n, sizeof(double));
    for(i = 0; i < dim_n; i++) {
      ALocalR[i] = ALocal((dim_n-2), i);
      VLocalR[i] = VLocal((dim_n-2), i);
    }
    
    for(i = P_ROWS-2; i > 0; i=i-2) {
      for(j = 0; j < dim_n; j++) {
        if(i == 2) {
          ALocal(i, j) = ALocal((i-1), j);
          VLocal(i, j) = VLocal((i-1), j);
        }
        else {
          ALocal(i, j) = ALocal((i-2), j);
          VLocal(i, j) = VLocal((i-2), j);
        }
      }
    }
    
    for(i = 1; i < P_ROWS-2; i=i+2) {
      for(j = 0; j < dim_n; j++) {
        ALocal(i, j) = ALocal((i+2), j);
        VLocal(i, j) = VLocal((i+2), j);
      }
    }    
    // assign rightmost-1 row to rightmost row
    for(i = 0; i < dim_n; i++) {
      ALocal((dim_n-1),i) = ALocalR[i];
      VLocal((dim_n-1),i) = VLocalR[i];
    }
  }
}


int main(int argc, char *argv[]) 
{
  
  int k, p;
  double start_time, end_time, elapsed_time;
  sumLocal = (double *) calloc(2, sizeof(double));
  // input parameter
  srand48(1);
  dim_n = atoi(argv[1]);
  max_sweeps = atoi(argv[2]);
  tol = pow(10, atoi(argv[3]));
  MPI_flag = atoi(argv[4]);
  N = dim_n*dim_n;
  threshold = tol * dim_n;

  // MPI init
  if(MPI_flag == 1) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);  
  }
  else {
    npes = 1;
  }
  
  P_ROWS = dim_n / npes;
  num_el = N / npes;
  sum_gap = log(dim_n)/log(2);
  if(sum_gap < 8) sum_gap = 8;
  
	start_time = MPI_Wtime();

  // Generate A V ALocal VLocal  
  generateALocalVLocal();
  for(iter_n = 0; (iter_n < max_sweeps) && (sum > threshold); iter_n++) {
    sumLocal[0] = 0; sumLocal[1] = 0;
    if(MPI_flag == 1) {
      // Create the row_type
      MPI_Type_contiguous(dim_n, MPI_DOUBLE, &row_type);
      MPI_Type_commit(&row_type);
      MPI_Type_size(row_type,&row_size);
    }
    // repeat individual Brenty-Luk moves 2*npes-1 times
    for (k=0;k<2*npes-1;k++)  {
      for(p=0; p<P_ROWS/2; p++)
        rotationComputation(ALocal+p*2*dim_n, VLocal+p*2*dim_n);  
      switchRows(); 
    }
    
    // sum sumLocal from all PEs and then update sum and threshold_tmp
    if(iter_n % sum_gap == 0){
      if(MPI_flag == 1){
        MPI_Allreduce(&sumLocal[0], &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&sumLocal[1], &threshold_tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      else {
        sum = sumLocal[0];
        threshold_tmp = sumLocal[1];
      }
      threshold = threshold_tmp * dim_n * tol;
      if(sum < threshold && rank == 0)
        printf("after iter %d times, sum = %.4e < threshold = %.4e\n", iter_n, sum, threshold);
    }
  }
  
  // get A
  if(MPI_flag == 1) {
    // gather ALocal and VLocal to get A and V
    MPI_Gather(ALocal, num_el, MPI_DOUBLE, A, num_el, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(VLocal, num_el, MPI_DOUBLE, V, num_el, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else {
    for(i = 0; i<dim_n; i++){ 
      for(j=0;j<dim_n;j++) {
          A(i,j) = ALocal(i,j);
          V(i,j) = VLocal(i,j);
      }
    }
  }
  
  // get V and Sigma
  if(rank == 0 || MPI_flag == 0) {

    // get the V^T from V
    for(i=0; i<dim_n; i++){ 
      for(j=i;j<dim_n;j++) {
          double tmp = V(i,j);
          V(i,j) = V(j,i);
          V(j,i) = tmp;
      }
    }  
    // get B = A^T*A 
    double* B = (double *) malloc(N*sizeof(double));
    double* Sigma = (double *) malloc(N*sizeof(double));
    for(i = 0; i < dim_n; i++) {
      for(j = 0; j < dim_n; j++) {
        B(i, j) = 0;
        for(k = 0; k < dim_n; k++) {
          B(i, j) += A(i, k) * A(j, k);
        }
      }
    }
    // get Sigma from B
    for (i=0;i<dim_n;i++){
 	    Sigma(i, i) = sqrt(B(i, i));
    }  
    
    // print Sigma and V
    /*  
    printf("rank %d Sigma after Brenty-Luk moves\n",rank);
    for (i=0;i<dim_n;i++){
      for (j=0;j<dim_n;j++){
    	  printf("%f ",Sigma(i, j));
      }
      printf("\n");
    }  

    printf("rank %d V after Brenty-Luk moves\n",rank);
    for (i=0;i<dim_n;i++){
      for (j=0;j<dim_n;j++){
    	  printf("%f ",V(i, j));
      }
      printf("\n");
    }
    */
    
    // get time
    end_time = MPI_Wtime();
    elapsed_time = end_time - start_time;
  	printf("dim_n = %d, npes = %d, time elapsed = %.4e sec\n", dim_n, npes, elapsed_time);    
  }
  
  // if use MPI funialize MPI
  if(MPI_flag == 1)
    MPI_Finalize();
  return 0;

  

}
