#! /bin/sh


# set CUDA and MPI environment path 
CUDA_HOME=/usr/local/cuda-10.2	# change this with your own configuration
MPI_HOME=/home/wyf/openmpi4		# change this with your own configuration

INC=-I$(CUDA_HOME)/include -I$(MPI_HOME)/include 
LIB=-L$(CUDA_HOME)/lib64 -L$(MPI_HOME)/lib

# set CUDA and MPI Dynamic link library
LINK= -lcudart -lcufft -lm -lmpi -lpthread -lrt -DMPICH_IGNORE_CXX_SEEK  -DMPICH_SKIP_MPICXX -lstdc++


all:
	@echo "***************************************************************************";
	@echo "**** Create a file folder for saving synthetic seismograms ****************";
	@if [ ! -d "./input/records_marmousi" ]; then \
		mkdir ./input/records_marmousi ; \
	fi
	@echo "The ./input/records_marmousi for saving Marmousi expQRTM records";

	@if [ ! -d "./input/records_crosswell" ]; then \
		mkdir ./input/records_crosswell ; \
	fi
	@echo "The ./input/records_crosswell for saving Crosswell expQRTM records"; 

	@echo "***************************************************************************";
	@echo "**** Create a file folder for saving output seismograms and images ********";
	@if [ ! -d "./output/output_marmousi" ]; then \
		mkdir ./output/output_marmousi ; \
	fi
	@echo "The ./output/output_marmousi for saving Marmousi expQRTM results";

	@if [ ! -d "./output/output_crosswell" ]; then \
		mkdir ./output/output_crosswell ; \
	fi
	@echo "The ./output/output_crosswell for saving Crosswell expQRTM results"; 


# CUDA and C++ source codes for Marmousi model
CFILES_marm = QRTM_marmousi.cpp CQRTM.cpp
CUFILES_marm = CUDAQRTM.cu
OBJECTS_marm = QRTM_marmousi.o CUDAQRTM.o 
EXECNAME_marm = marmousi.out


marmousi:
	mpicc -w -c $(CFILES_marm) $(INC) $(LIB) $(LINK) 
	nvcc -w -c $(CUFILES_marm) $(INC) $(LIB) $(LINK) 
	mpicc -o $(EXECNAME_marm) $(OBJECTS_marm) $(INC) $(LIB) $(LINK) 

	rm -f *.o 
	nohup ./$(EXECNAME_marm) &


# CUDA and C++ source codes for Crosswell model
CFILES_cw = QRTM_crosswell.cpp CQRTM.cpp
CUFILES_cw = CUDAQRTM.cu
OBJECTS_cw = QRTM_crosswell.o CUDAQRTM.o 
EXECNAME_cw = crosswell.out


crosswell:
	mpicc -w -c $(CFILES_cw) $(INC) $(LIB) $(LINK) 
	nvcc -w -c $(CUFILES_cw) $(INC) $(LIB) $(LINK) 
	mpicc -o $(EXECNAME_cw) $(OBJECTS_cw) $(INC) $(LIB) $(LINK) 

	rm -f *.o 
	nohup ./$(EXECNAME_cw) &