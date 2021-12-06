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
├── input
│   ├── para_crosswell
│   ├── para_marmousi
│   ├── records_crosswell
│   └── records_marmousi
├── output
│   ├── output_crosswell
│   └── output_marmousi
├── examples
│   ├── plot_crosswell
│   └── plot_marmousi
├── Myfunctions.h
├── CQRTM.cpp
├── CUDAQRTM.cu
├── QRTM_crosswell.cpp
├── QRTM_marmousi.cpp
├── Makefile
└── README.md
```

### The input parameters of the package

`/input/para_*`: accurate velocity and $Q$ model for Marmousi and crosswell examples:
    - `/input/para_marmousi/acc_vp.dat`: quasi-Marmousi velocity model;
    - `/input/para_marmousi/acc_Qp.dat`: quasi-Marmousi $Q$ model;
    - `/input/para_crosswell/acc_vp.dat`: synthetic crosswell $Q$ model;
    - `/input/para_crosswell/acc_Qp.dat`: synthetic quasi-Marmousi $Q$ model;
    - `/input/para_crosswell/sx.txt`: the x location of sources in grid;
    - `/input/para_crosswell/sz.txt`: the z location of sources in grid;
    - `/input/para_crosswell/rx.txt`: the x location of receivers in grid;
    - `/input/para_crosswell/rz.txt`: the z location of receivers in grid.



### The output results of the package    

`output/output_*`: generated results such as seismograms, images of each shot and final stacked images. What we are most interested in is final images, wich includes:
    - `output/output_*/Final_image_cor_type0.dat`: image from acoustic RTM;
    - `output/output_*/Final_image_cor_type1.dat`: image from viscoacoustic RTM without compensation;
    - `output/output_*/Final_image_cor_type2.dat`: image from $Q$-RTM without stabilization;
    - `output/output_*/Final_image_cor_type3.dat`: image from $Q$-RTM using explicit stabilization (alpha=1);
    - `output/output_*/Final_image_cor_type4.dat`: image from $Q$-RTM using explicit stabilization (alpha=2);
    - `output/output_*/Final_image_cor_type5.dat`: image from $Q$-RTM using explicit stabilization (alpha=8);


### The examples for plotting the final results    

`examples/plot_*`: scripts for plotting final results, which includes:
    - `examples/plot_*/data`: a folder for final results copied from `./output/output_*`;
    - `examples/plot_*/madaimage/SConstruct`: plot images, velocity and $Q$ models;
    - `examples/plot_*/matlabtrace/martrace`: plot extracted trace from final migrated images for comparison.

### The C and CUDA files    

-   `Myfunctions.h`: header file for the statement of extern "C" functions defined in CUDAQRTM.cu and C functions defined in CQRTM.cpp;
-   `CQRTM.cpp`: C functions file;
-   `CUDAQRTM.cu`: CUDA functions file;
-   `QRTM_*.cpp`: C main function file;
-   `Makefile`: excution script.



## Important flags and parameters in QRTM_*.cpp

- `Geometry`: Geometry=0 for surface seismic survey and Geometry=1 for crosswell seismic survey;
- `RTMtype`: you can change this flag to generate different migrated images:
    
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

- explicit satbiliztion parameters for different `RTMtype` configurations:
    
``` c
	float alphaorder = 2;	// the default stabilization order;
	float kref = 0.32;	// the reference wavenumber for stabilization;
	float scaling = 0;	// the scaling factor see equation (23) in our paper;
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

# set CUDA and MPI environment path 
CUDA_HOME=/usr/local/cuda-10.2	# change this with your own configuration
MPI_HOME=/home/wyf/openmpi4		# change this with your own configuration

INC=-I$(CUDA_HOME)/include -I$(MPI_HOME)/include 
LIB=-L$(CUDA_HOME)/lib64 -L$(MPI_HOME)/lib
```

- Step 2: Run the `Makefile` by the command line: `make` for generating records folder for saving synthetic seismograms;
- Step 3: Run `make marmousi` for surface seismic imaging on Marmousi model, or run `make crosswell` for crosswell imaging;
- Step 4: View generated files in the folder `./ouput/output_*`;
- Step 5: Plot final results by run `/examples/plot_*/madaimage/SConstruct` and `/examples/plot_*/matlabtrace/martrace.m`.


## Contact me

I am Yufeng Wang, an associate research professor from China University of Geosciences, Wuhan. If you have any question about this coda package, please feel free to contact me by [Email:wangyufeng@cug.edu.cn](wangyufeng@cug.edu.cn).

## Copyright

`expQRTM` is a CUDA-based code package that implements Q-compensated reverse time migration with explicit stabiization in the time-space domain. The Marmousi model used in this package is available for download from `Madagascar` MainPage <http://www.ahay.org/data/marm2>.

Copyright (C) 2021  China University of Geosciences, Wuhan (Yufeng Wang)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTI CULAR PURPOSE.  See theGNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
