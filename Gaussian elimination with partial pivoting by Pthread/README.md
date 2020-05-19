# 1. code  
The code is in the file of gaussian.c.  
# 2. logic  
## (1) Pipeline  
 
![flow](https://github.com/ZhixinLai/Parallel-Computing/blob/master/Gaussian%20elimination%20with%20partial%20pivoting%20by%20Pthread/flow.png)    

<p align="center">Fig1 Logic flow chart</p>  
## (2) sub-processes (how to partition the data and how they are synchronized)  
### •	1.1 Find the max value in A(i:n, i) and is m   
Using p threads to find the max value in A(i:n, i). Each thread compares some values with the global max value and update the global max variable. If (n-i)%p == 0, each thread get (n-i)/p works, otherwise most threads get ((n-i)/p+1) works and one gets the left works. In order to avoid the race condition, I use the locking method with pthread_mutex_lock function. After all the thread finish their work, I use pthread_join function to get all the thread synchronized. The GetMaxValue in the appendix shows more details.  
### •	1.3 Update A(i+1:n,i+1:n) and b(i+1:n)  
Using p threads to Update A(i+1:n,i+1:n) and b(i+1:n) and get a up triangle matrix A. Each thread is to calculate the (n-i-1)/p lines. If (n-i-1) % p == 0, each thread get (n-i-1)/p works, otherwise most threads get ((n-i-1)/p+1) works and one gets the left works. Because the values of line i are fixed and each thread just adjust local data, so the work among different threads are independent. After all the thread finish their work, I use pthread_join function to get all the thread synchronized. There is not need for locking function. The GetUpTriangle in the appendix shows more details.  
### •	2 Update A(0:i,i) and b(0:i)  
Using p threads to Update A(0:i,i) and b(0:i) and get a diagonal matrix A. Each thread is to calculate the i/p lines. If i%p == 0, each thread get i/p works, otherwise most threads get (i/p+1) works and one gets the left works. Each thread works independently, so there is not need for locking function. After all the thread finish their work, I use pthread_join function to get all the thread synchronized. The GetDiagonal in the appendix shows more details.  
### •	3 Get the final result  
Using p threads to calculate the result. Each thread is to calculate the dim_n/p lines and each line is just do a division B[i]/A[i][i]. If n%p == 0, each thread get n/p works, otherwise most threads get (n/p+1) works and one gets the left works. Each thread works independently, so there is not need for locking function. After all the thread finish their work, I use pthread_join function to get all the thread synchronized. The GetResultin the appendix shows more details.  
### •	Others  
Other parts don’t include any parallel computing method.  

# 3. Benchmark  
First, I test the cores of the server first by using the command below.  
cat /proc/cpuinfo| grep "physical id"| sort| uniq| wc -l  
cat /proc/cpuinfo| grep "cpu cores"| uniq  
The first get the CPU number (10) and the second get the core number of each CPU(1). Therefore, the total core number is 10 * 1 = 10.  
According to the requirement, the start n = 512 and increase it by a factor of 2 till n = 4096 and start p = 1 and increase to twice the number of cores. Therefore, the dim_n = 512, 1024, 2048, 4096, and p = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20.  

# 4. Result  
The results are shown in the Fig2, the (a)(b)(c)(d) shows the results with dimension of A being 512, 1024, 2048, 4096 respectively. For each graph, the x is the number of threads from 1 to 20 and the y means the total time of solving the equation.  
     

![flow](https://github.com/ZhixinLai/Parallel-Computing/blob/master/Gaussian%20elimination%20with%20partial%20pivoting%20by%20Pthread/result.png)    

<p align="center">Fig2 result</p>  

# 5. Discussion  
### •	Influence of the thread number  
As we can see from the Fig2, when the thread number grows, the consuming time gets lower first and grows then. The reason for time decreasing with thread growing is that each thread can work together and the work load of each thread is less than the total work load. It is easy to understand. However, when the thread number increases further, the consuming time gets higher, I guess the reason is that the process of assigning task to different threads, creating threads, joining threads, locking and other thread related manipulations consumes time. We can call the time cost as communication cost. When the thread number grows, the communication cost gets higher. Once the communication overhead compensates more than the time saved by multi-thread computing, the total time will grow with p grows.   
### •	Influence of the matrix dimension  
First conclusion: when the p set the same value, the run time grows with the n growing. It is easy to understand. If the n grows, the floating computing number will grow.  
Second conclusion: when n grows, the local minimum point moves in the right direction and the growing trend after the local minimum point is lower. The p values of local minimum points are 3 4 4 8 for matrix dimension 512 1024 2048 4096 respectively. As we discussion in the former part, communication overhead time is Tw and float computing time is Tc.  The ratio of Tw to Tc is high when the matrix dimension is low, while the ratio is low when the matrix dimension is high. Therefore, when n = 4096, the Tw takes a small part of the total time(Tw + Tc) and thus the increasing trend is low when Tw grows.   
