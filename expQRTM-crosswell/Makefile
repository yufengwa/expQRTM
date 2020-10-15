#! /bin/sh

# compiler

# set CUDA and MPI environment path 
CUDA_HOME=/usr/local/cuda-10.0
MPI_HOME=/home/wyf/mpich3

MPICC = $(MPI_HOME)/bin/mpicc

INC=-I$(CUDA_HOME)/include -I$(MPI_HOME)/include 
LIB=-L$(CUDA_HOME)/lib64 -L$(MPI_HOME)/lib -L/usr/lib/x86_64-linux-gnu

# set CUDA and MPI Dynamic link library
LINK= -lcudart -lcufft -lm -lmpich -lpthread -lrt -DMPICH_IGNORE_CXX_SEEK  -DMPICH_SKIP_MPICXX -lstdc++

# CUDA and C++ source codes
CFILES = QRTM.cpp
CUFILES = CUDAQRTM.cu
OBJECTS = QRTM.o CUDAQRTM.o 
EXECNAME = a.out

all:
	mpicc -w -c $(CFILES) $(INC) $(LIB) $(LINK) 
	nvcc -w -c $(CUFILES) $(INC) $(LIB) $(LINK) 
	mpicc -o $(EXECNAME) $(OBJECTS) $(INC) $(LIB) $(LINK) 

	rm -f *.o 
	nohup ./a.out&
