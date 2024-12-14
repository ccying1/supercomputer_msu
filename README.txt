#BSUB -n 1
#BSUB -W 00:15
#BSUB -o "task1_4040_16.%J.out"
#BSUB -e "task1_4040_16.%J.err"
#BSUB -R "span[hosts=1]"
OMP_NUM_THREADS=16 ./task1_4040

#BSUB -n 1
#BSUB -W 00:15
#BSUB -o "task1_4040_4.%J.out"
#BSUB -e "task1_4040_4.%J.err"
#BSUB -R "span[hosts=1]"
OMP_NUM_THREADS=4 ./task1_4040

#BSUB -n 1
#BSUB -W 00:15
#BSUB -o "task1_4040_1.%J.out"
#BSUB -e "task1_4040_1.%J.err"
#BSUB -R "span[hosts=1]"
OMP_NUM_THREADS=1 ./task1_4040



#!/bin/bash
#BSUB -J task2_2mpi          # 作业名称
#BSUB -o task2_2mpi_output.%J    # 输出文件（%J 表示作业ID）
#BSUB -e task2_2mpi_error.%J     # 错误文件
#BSUB -n 4                # 申请的核心数
#BSUB -R "span[ptile=2]"  # 每个节点的核心数
#BSUB -q normal           # 指定队列名，根据系统配置填写
#BSUB -W 00:30            # 预计运行时间 (格式: HH:MM)

# 加载必要的模块 (根据系统实际环境调整)
module load OpenMPI/4.0.0
mpicxx -std=c++11 -o task2mpi_2_4040 task2mpi.cpp -lm
# 执行 MPI 程序
mpirun -np 2 ./task2mpi_2_4040


#!/bin/bash
#BSUB -J task2_4mpi          # 作业名称
#BSUB -o task2_4mpi_output.%J    # 输出文件（%J 表示作业ID）
#BSUB -e task2_4mpi_error.%J     # 错误文件
#BSUB -n 4                # 申请的核心数
#BSUB -R "span[ptile=2]"  # 每个节点的核心数
#BSUB -q normal           # 指定队列名，根据系统配置填写
#BSUB -W 00:30            # 预计运行时间 (格式: HH:MM)

# 加载必要的模块 (根据系统实际环境调整)
module load OpenMPI/4.0.0
mpicxx -std=c++11 -o task2mpi_4_4040 task2mpi.cpp -lm
# 执行 MPI 程序
mpirun -np 4 ./


#!/bin/bash
#BSUB -J task2_1mpi          # 作业名称
#BSUB -o task2_1mpi_output.%J    # 输出文件（%J 表示作业ID）
#BSUB -e task2_1mpi_error.%J     # 错误文件
#BSUB -n 4                # 申请的核心数
#BSUB -R "span[ptile=2]"  # 每个节点的核心数
#BSUB -q normal           # 指定队列名，根据系统配置填写
#BSUB -W 00:30            # 预计运行时间 (格式: HH:MM)

# 加载必要的模块 (根据系统实际环境调整)
module load OpenMPI/4.0.0
mpicxx -std=c++11 -o task2mpi_1_4040 task2mpi.cpp -lm
# 执行 MPI 程序
mpirun -np 1 ./task2mpi_1_4040




#!/bin/bash
#BSUB -J task3openmp_mpi24        # 作业名称
#BSUB -o task3openmp_mpi24_output.%J  # 标准输出文件
#BSUB -e task3openmp_mpi24_error.%J   # 错误输出文件
#BSUB -q short                # 提交到的队列名称
#BSUB -n 4                     # 请求的总核心数（MPI 进程数）
#BSUB -R "span[ptile=4]"       # 每个节点最多运行 4 个进程
#BSUB -x                       # 排他使用节点

# 加载所需的模块
module load OpenMPI/4.0.0      # 替换为实际的 OpenMPI 版本
module load GCC/9.3.0          # 替换为实际的 GCC 版本
# 设置 OpenMP 线程数
export OMP_NUM_THREADS=4       # 每个 MPI 进程的 OpenMP 线程数
# 编译代码（如果未提前编译）
mpicxx -fopenmp -std=c++11 -o task3openmp_mpi8090_2_4 task3openmp_mpi.cpp -lm
# 运行 MPI + OpenMP 程序
mpirun -np 2 ./task3openmp_mpi4040_2_4




#!/bin/bash
#BSUB -J task3openmp_mpi14        # 作业名称
#BSUB -o task3openmp_mpi14_output.%J  # 标准输出文件
#BSUB -e task3openmp_mpi14_error.%J   # 错误输出文件
#BSUB -q normal                # 提交到的队列名称
#BSUB -n 4                     # 请求的总核心数（MPI 进程数）
#BSUB -R "span[ptile=4]"       # 每个节点最多运行 4 个进程
#BSUB -x                       # 排他使用节点

# 加载所需的模块
module load OpenMPI/4.0.0      # 替换为实际的 OpenMPI 版本
module load GCC/9.3.0          # 替换为实际的 GCC 版本

# 设置 OpenMP 线程数
export OMP_NUM_THREADS=4       # 每个 MPI 进程的 OpenMP 线程数

# 编译代码（如果未提前编译）
mpicxx -fopenmp -std=c++11 -o task3openmp_mpi160180_1_4 task3openmp_mpi4040.cpp -lm

# 运行 MPI + OpenMP 程序
mpirun -np 1 ./task3openmp_mpi4040_1_4