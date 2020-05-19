Two One sided Jacobi for computing the EVD
==========================================

# 1. code

The code is in the file of jacobi.c, jacobi.sub, jacobi.sh.

# 2. logic

## (1) Pipeline

![flow](https://github.com/ZhixinLai/Parallel-Computing/blob/master/One%20sided%20Jacobi%20for%20computing%20the%20EVD%20with%20MPI/flow.png)    

<p align="center">Fig1 Logic flow chart</p>  

## (2) sub-processes (how to partition the data and how they are synchronized)  
### •	rotationComputation part  
I divide dim_n rows data into npes tasks. Within a task, I call the rotationComputation for P_ROW/2 times and each time just computation 2 rows in neighbor. Rows in different tasks don’t influence each other, therefore, it can be computed parallelly.  
### •	switchRows part  
In this part, rows should be exchanged among different tasks and within a task. As for exchange among different tasks, I use MPI_Send, MPI_Recv, and MPI_Barrier. More details can be seen in the code file. As for exchange within a task, I just divide them in to odd and even. Besides, the two situations should pay attention to special cases of rank = 0 and rank = npes – 1.  

# 3. Benchmark  
According to the requirement, the start n = 128 and increase it by a factor of 2 till n = 1024 and start task number(npes) = 1 and increase to twice the number of cores. Therefore, the dim_n = 128, 256, 512, 1024 and npes = 1 2 4 8 16 32 64.  
Besides, in order to compare the achieved speed-ups versus a sequential code (not an MPI code), I extend another npes = 0. npes = 0 means that no MPI is used in the function, of which the result represents the performance of a sequential code.  

# 4. Result  
The results are shown in the Fig3, it shows the results with dimension dim_n being 128, 256, 512, 1024 respectively. For each graph, the x is the number of task number npes from 0 to 64 and the y means the total run time.  
![flow](https://github.com/ZhixinLai/Parallel-Computing/blob/master/One%20sided%20Jacobi%20for%20computing%20the%20EVD%20with%20MPI/result.png)  
<p align="center">Fig2 Results</p>  

# 5. Discussion  
### •	Comparison between sequential computing and parallel computing  
As we can see from the result show before. npes = 0 means the code runs in sequential way without MPI method. P = 1 means that we use MPI with npes = 1. No matter how many dim_n is, the computing time of npes = 1 is greater than npes = 0, which means although they both use 1 thread, the extra manipulations for MPI cost time. After comparing the npes = 0 and npes > 1, the efficiency of parallel computing is better than sequential computing in most situations especially when dim_n is large(256 512 1024). When dim_n is small(128), the efficiency of parallel computing is worse than sequential computing.  
### •	Influence of the task number (npes > 0)  
As we can see from the Fig2, when the npes grows, the consuming time gets lower first and grows then. The reason for time decreasing with npes growing is that each process can work together and the work load of each process is less than the total work load. It is easy to understand. However, when the npes increases further, the consuming time gets higher (not obvious but does exist especially when dim_n is low), the reason is that the process of assigning task to different processes, data transmission among different processes, and barrier manipulations consume time. We call the time cost as extra cost. When the npes grows, the extra cost gets higher. Once the extra cost compensates more than the time saved by multi-processes computing, the total time will grow with npes grows. Besides, when dim_n get higher, the increase trend gets lower, the reason is that the ratio of extra cost to computing cost get lower when dim_n get larger.   
### •	Influence of the matrix dimension  
First conclusion: when the npes set the same value, the run time grows with the dim_n growing. It is easy to understand. If the dim_n grows, the floating computing number will grow.  
Second conclusion: when dim_n grows, the local minimum point moves in the right direction and the growing trend after the local minimum point is lower. As we discussion in the former part, extra cost is Tw and float computing time is Tc.  The ratio of Tw to Tc is high when the matrix dimension is low, while the ratio is low when the matrix dimension is high.   

