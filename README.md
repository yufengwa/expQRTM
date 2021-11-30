---
title: expQRTM Manual
date: 2021-11
author: Yufeng Wang
mathjax: true
---


## Overview of `expQRTM` package

`expQRTM` is a CUDA-based code package that implements $Q$-compensated reverse time migration with explicit stabiization in the time-space domain. This package is provided for verifying the feasibility and stability of the proposed explicit stabilization scheme. We provide two package versions: `QRTM_marmousi.cpp` for $Q$-RTM on a surface survey with Marmousi model; `QRTM_crosswell.cpp` for $Q$-RTM on a crosswell survey. The package has been well test on Linux OS (Ubuntu 20.04 TLS) with `MPI` and `CUDA` available.

## The architecture of `expQRTM` package 

``` bash
├── CUDAQRTM.cu
├── examples
│   ├── plot_crosswell
│   └── plot_marmousi
├── input
│   ├── para_crosswell
│   └── para_marmousi
├── Makefile
├── Myfunctions.h
├── output
│   ├── output_crosswell
│   └── output_marmousi
├── QRTM_crosswell.cpp
├── QRTM_marmousi.cpp
└── README.md
```


-   `/input/para_*`: accurate velocity and $Q$ model for $Q$-RTM:
    - `acc_vp.dat`: quasi-Marmousi velocity model;
    - `acc_Qp.dat`: quasi-Marmousi $Q$ model.
    
-   `output/output_*`: generated results such as seismograms, images of each shot and final stacked images. What we are most interested in is final images, wich includes:
    - `Final_image_cor_type0.dat`: image from acoustic RTM;
    - `Final_image_cor_type1.dat`: image from viscoacoustic RTM without compensation;
    - `Final_image_cor_type2.dat`: image from $Q$-RTM without stabilization;
    - `Final_image_cor_type3.dat`: image from $Q$-RTM using explicit stabilization (alpha=1);
    - `Final_image_cor_type4.dat`: image from $Q$-RTM using explicit stabilization (alpha=2);
    - `Final_image_cor_type5.dat`: image from $Q$-RTM using explicit stabilization (alpha=8);
    
-   `examples/plot_*`: scripts for plotting figures, which includes:
    - `/data`: a folder for final results copied from `./output`;
    - `/madaimage/SConstruct`: plot images, velocity and $Q$ models;
    - `/matlabtrace/martrace`: plot extracted trace from final migrated images for comparison.
    
-   `Myfunctions.h`: header file;
-   `CUDAQRTM.cu`: cuda code file;


-   `QRTM_*.cpp`: c++ code file, there are serveal important flags and parameters to control performance of $Q$-RTM, which includes:

    - `RTMtype`: you can change this flag to generate different migrated images.
    
``` c
	int Geometry=0;		// Geometry=0 for surface seismic survey
					// Geometry=1 for crosswell seismic survey
						
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

- CUDA environment (cuda-10.2 for example);
- MPI environment (openmpi4 for example);
- Matlab;
- Madagascar.


## How to run this package

- Step 1: Confirm the environment in `Makefile`, and replace the folder path with your own enviroment path; 

``` bash
#! /bin/sh

# compiler

# set CUDA and MPI environment path 
CUDA_HOME=/usr/local/cuda-10.2
MPI_HOME=/home/wyf/openmpi4

INC=-I$(CUDA_HOME)/include -I$(MPI_HOME)/include 
LIB=-L$(CUDA_HOME)/lib64 -L$(MPI_HOME)/lib

# set CUDA and MPI Dynamic link library
LINK= -lcudart -lcufft -lm -lmpi -lpthread -lrt -DMPICH_IGNORE_CXX_SEEK  -DMPICH_SKIP_MPICXX -lstdc++

```
- Step 2: Run the `Makefile` by the command line: `make` for generating records folder for synthetic seismograms;
- Step 3: Run `make marmousi` for surface seismic imaging on Marmousi model, or run `make crosswell` for crosswell imaging;
- Step 4: View generated files in the folder `./ouput/output_*`;
- Step 5: Plot figures by run `/examples/plot_*/madaimage/SConstruct` and `/examples/plot_*/matlabtrace/martrace.m`.


## Contact me

I am Yufeng Wang, an associate researcher from China University of Geosciences, Wuhan. If you have any question about this coda package, please feel free to contact me by [Email:wangyufeng@cug.edu.cn](wangyufeng@cug.edu.cn).

## Copyright

`expQRTM` is a CUDA-based code package that implements Q-compensated reverse time migration with explicit stabiization in the time-space domain. The Marmousi model used in this package is available for download from `Madagascar` MainPage <http://www.ahay.org/data/marm2>

Copyright (C) 2021  China University of Geosciences, Wuhan (Yufeng Wang)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTI CULAR PURPOSE.  See theGNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
