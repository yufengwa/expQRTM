---
title: exp$Q$RTM Manual
date: 2020-10
author: Yufeng Wang
mathjax: true
---


## Overview of `expQRTM` package

`expQRTM` is a CUDA-based code package that implements $Q$-compensated reverse time migration with explicit stabiization in the time-space domain. This package is provided for verifying the feasibility and stability of the proposed explicit stabilization scheme. We provide two package versions: `expQRTM-marm` for $Q$-RTM on a surface survey with Marmousi model; `expQRTM-crosswell` for $Q$-RTM on a crosswell survey. The package has been well test on Linux OS (Ubuntu 18.04 TLS) with `MPI` and `CUDA` available.

## The architecture of `expQRTM` package 

-   `input`: accurate velocity and $Q$ model for $Q$-RTM:
    - `acc_vp.dat`: quasi-Marmousi velocity model;
    - `acc_Qp.dat`: quasi-Marmousi $Q$ model;
    - `./acorecords`: a folder for acoustic records after removing first arrivals;
    - `./attrecords`: a folder for attenuated records after removing first arrivals;
-   `output`: generated results such as seismograms, images of each shot and final stacked images. What we are most interested in is final images, wich includes:
    - `Final_image_cor_type0.dat`: image from acoustic RTM;
    - `Final_image_cor_type1.dat`: image from viscoacoustic RTM without compensation;
    - `Final_image_cor_type2.dat`: image from $Q$-RTM without stabilization;
    - `Final_image_cor_type3.dat`: image from $Q$-RTM using explicit stabilization (alpha=1);
    - `Final_image_cor_type4.dat`: image from $Q$-RTM using explicit stabilization (alpha=2);
    - `Final_image_cor_type5.dat`: image from $Q$-RTM using explicit stabilization (alpha=8);
-   `plot`: scripts for plotting figures, which includes:
    - `./data`: a folder for final results copied from `./output`;
    - `/madaimage/SConstruct`: plot images, velocity and $Q$ models;
    - `/matlabtrace/martrace`: plot extracted trace from final migrated images for comparison.
-   `Myfunctions.h`: header file;
-   `CUDAQRTM.cu`: cuda code file;
-   `QRTM.cpp`: c++ code file, there are serveal important flags and parameters to control performance of $Q$-RTM, which includes:
    - `RTMtype`: you can change this flag to generate different migrated images.
``` c
	int RTMtype=0;	// RTMtype=0 for acoustic RTM
			// RTMtype=1 for viscoacoustic RTM without compensation
			// RTMtype=2 for QRTM withou stabilization
			// RTMtype=3 for QRTM using explicit stabilization scheme (alpha=1)
			// RTMtype=4 for QRTM using explicit stabilization scheme (alpha=2)
			// RTMtype=5 for QRTM using explicit stabilization scheme (alpha=8)
```
    - parameters for explicit satbiliztion scheme:
``` c
	float alphaorder = 2;
	float kref = 0.32;
	float scaling = 0;
	if(RTMtype==2)
	{
		scaling = 0;
	}
	if(RTMtype==3)
	{
		alphaorder = 1;
		scaling = 1;
	}
	if(RTMtype==4)
	{
		alphaorder = 2;
		scaling = 1;
	}
	if(RTMtype==5)
	{
		alphaorder = 8;
		scaling =1;
	}
	float sigmafactor = scaling*powf(kref, 2*averGamma+1-alphaorder);
```
-   `Makefile`: excution script.


## Prerequisites

`expQRTM` package is developed under `Linux` system, which should be equipped with the following environments:

- CUDA environment (for example, `-I/usr/local/cuda-10.0/include` `-L/usr/local/cuda-10.0/lib64`);
- MPI environment (for example, `-I/home/wyf/intel/impi/5.0.1.035/intel64/include` `-L/home/wyf/intel/impi/5.0.1.035/intel64/lib/`);
- matlab;
- madagascar.


## How to run this package

- Step 1: Confirm the environment in `Makefile`, and replace the folder path with your own enviroment path; 

``` bash
#! /bin/sh
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
# Execution
all:
	mpicc -w -c $(CFILES) $(INC) $(LIB) $(LINK) 
	nvcc -w -c $(CUFILES) $(INC) $(LIB) $(LINK) 
	mpicc -o $(EXECNAME) $(OBJECTS) $(INC) $(LIB) $(LINK) 
	rm -f *.o 
	nohup ./a.out&
```
- Step 2: Run the `Makefile` by the command line: `make`;
- Step 3: View generated files in the folder `./ouput`;
- Step 4: Plot figures by run `/plot/madagascar/SConstruct` and `/plot/matlab/martrace.m`.


## Contact me

I am Yufeng Wang, an associate researcher from China University of Geosciences, Wuhan. If you have any question about this coda package, please feel free to contact me by [Email:wangyufeng@cug.edu.cn](wangyufeng@cug.edu.cn).

## Copyright

`expQRTM` is a CUDA-based code package that implements Q-compensated reverse time migration with explicit stabiization in the time-space domain. The Marmousi model used in this package is available for download from `Madagascar` MainPage <http://www.ahay.org/data/marm2>

Copyright (C) 2020  China University of Geosciences, Wuhan (Yufeng Wang)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTI CULAR PURPOSE.  See theGNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
