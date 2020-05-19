# 1. code  
The code is in the file of gaussian.c.
# 2. logic  
     (1) Pipeline
 
Fig1 Logic flow chart
## (2) sub-processes (how to partition the data and how they are synchronized)
### •	1.1 Find the max value in A(i:n, i) and is m
If p == 0, the getMaxwithoutOpenMP() will be called, which means that I do not use OpenMP to do parallel computing, otherwise the getMaxwithOpenMP() will be called. Using p threads to find the max value in A(i:n, i). Each thread compares some values with the global max value and update the global max variable. In order to avoid the race condition, I use the #pragma omp critical(name). I use the dynamic type for schedule, therefor the task among each thread will be even distributed. The getMaxwithOpenMP() function shows more details how I realize the parallel computing. 
### •	1.2 Swap A(:, i) with A(:, m); Swap B(:, i) with B(:, m)
The swapwithMax() function will be called in the step. I do not use parallel computing method. Just do it sequentially.
### •	1.3 Update A(i+1:n,i+1:n) and B(i+1:n, :)
If p == 0, the getUpTrianglewithoutOpenMP() will be called, which means that I do not use OpenMP to do parallel computing, otherwise the getUpTrianglewithOpenMP() will be called. Using p threads to update A(i+1:n,i+1:n) and B(i+1:n, :) and get a up triangle matrix A. Because the values of row i are fixed and each thread just adjust local data, so the work among different threads are independent. And thus I do not need for locking function. I use the dynamic type for schedule, therefor the task among each thread will be even distributed. The getUpTrianglewithOpenMP() function shows more details how I realize the parallel computing. 
### •	2 backsubbstitution
In this step, I divide it into two steps: Get the result X(i, :) = B(i, :) / A(i, i) and Update the B(1: i-1, :). If p == 0, the backSubbstitutionwithoutOpenMP() will be called, which means that I do not use OpenMP to do parallel computing, otherwise the backSubbstitutionwithOpenMP() will be called. Each thread works independently, so there is not need for locking function. I use the dynamic type for schedule, therefor the task among each thread will be even distributed. The backSubbstitutionwithOpenMP() shows more details how I realize the parallel computing. 

# 3. Benchmark
First, I test the cores of the server first by using the command below.
cat /proc/cpuinfo| grep "physical id"| sort| uniq| wc -l
cat /proc/cpuinfo| grep "cpu cores"| uniq
The first get the CPU number (10) and the second get the core number of each CPU(1). Therefore, the total core number is 10 * 1 = 10.
According to the requirement, the start n = 512 and increase it by a factor of 2 till n = 4096 and start p = 1 and increase to twice the number of cores. Therefore, the dim_n = 512, 1024, 2048, 4096, and p = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20.
Besides, in order to compare the achieved speed-ups versus a sequential code (not an openMP code with a single thread), I extend the input from 1-20 to 0-20. P = 0 means that no openMP is used in the function, of which the result represents the performance of a sequential code.

# 4. Result
The results are shown in the Fig2, the (a)(b)(c)(d) shows the results with dimension of A being 512, 1024, 2048, 4096 respectively. For each graph, the x is the number of threads from 0 to 20 and the y means the total time of solving the equation.
  

Fig2 Results
# 5. Discussion
### •	Comparison between sequential computing and parallel computing
As we can see from the result show before. P = 0 means the code runs in sequential way without OpenMP method. P = 1 means that the thread we use is 1 but use OpenMP. No matter how many dim_n is, the computing time of p = 1 is greater than p = 0, which means although they both use 1 thread, the extra manipulations for OpenMP cost time. After comparing the p = 0 and p > 1, the efficiency of parallel computing is better than sequential computing in most situations. There are some abnormal cases that computing time sharp increases, like the dim_n = 512 and p = 3.
### •	Influence of the thread number
As we can see from the Fig2, when the thread number grows, the consuming time gets lower first and grows then. The reason for time decreasing with thread growing is that each thread can work together and the work load of each thread is less than the total work load. It is easy to understand. However, when the thread number increases further, the consuming time gets higher, I guess the reason is that the process of assigning task to different threads, creating threads, locking and other thread related manipulations consumes time. We can call the time cost as extra cost. When the thread number grows, the extra cost gets higher. Once the communication overhead compensates more than the time saved by multi-thread computing, the total time will grow with p grows. Besides, when dim_n get higher, the increase trend gets lower, the reason is that the ratio of extra cost to computing cost get lower when dim_n get larger. 
### •	Influence of the matrix dimension
First conclusion: when the p set the same value, the run time grows with the n growing. It is easy to understand. If the n grows, the floating computing number will grow.
Second conclusion: when n grows, the local minimum point moves in the right direction and the growing trend after the local minimum point is lower. As we discussion in the former part, extra cost is Tw and float computing time is Tc.  The ratio of Tw to Tc is high when the matrix dimension is low, while the ratio is low when the matrix dimension is high. 
### •	Compare the result of P_thread and OpenMP
This part is not definitely true and is what I infer from the results of the experiments. I hope professor and TA could give some advice. From the result of the previous experiments. The curve when use OpenMP is less smooth than that when use pthread. The reason I guess is that, when I use pthread, I distribute the task evenly. When I use OpenMP, I use dynamic method to distribute task among different threads. I think it is a better way to distribute task evenly than static way, but I guess it can not still distribute it evenly, which may cost time of waiting each thread to join together. Do processor or TA have any ideas why there are more abnormal points when using OpenMP than using Pthread, like dim_n = 512 and p = 3?
