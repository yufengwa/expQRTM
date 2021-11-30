/* 2D GPU+MPI-based QRTM with explicit stabilization in the time-space domain for surface seismic imaging. */

/*
	Copyright (C) 2021  XXX (will be updated after acceptance)	    

	Acknowledgement:

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/


#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#define PI 3.1415926
using namespace std;
#include "Myfunctions.h"

int main(int argc, char* argv[])             
{
	//=========================================================
	//  MPI Indix
	//  =======================================================

	int myid,numprocs,namelen;	
	MPI_Comm comm=MPI_COMM_WORLD;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(comm,&myid);
	MPI_Comm_size(comm,&numprocs);
	MPI_Get_processor_name(processor_name,&namelen);
	if(myid==0)
	{
		printf("=========================================================\n");
		printf("=====================  expQRTM for Marmousi  ============\n");
		printf("=========================================================\n");
		printf("Number of MPI thread is %d\n",numprocs);
	}
		

	//=========================================================
	//  Parameters of the time of the system...
	//  =======================================================

	clock_t start, end;
	float runtime=0.0;


	//=========================================================
	//  Flags (in this package we set flags=0 as defult)
	//  =======================================================
	int Geometry=0;		// Geometry=0 for surface seismic survey
						// Geometry=1 for crosswell seismic survey

	int RTMtype=2;		// RTMtype=0 for acoustic RTM
						// RTMtype=1 for viscoacoustic RTM without compensation
						// RTMtype=2 for QRTM withou stabilization
						// RTMtype=3 for QRTM using explicit stabilization scheme (alpha=1)
						// RTMtype=4 for QRTM using explicit stabilization scheme (alpha=2)
						// RTMtype=5 for QRTM using explicit stabilization scheme (alpha=8)

	int Save_Not=0;		// Save_Not=1 for save forward, reconstruction and backward wavefileds
						// Save_Not=0 for don't save forward, reconstruction and backward wavefileds

	int Sto_Rec=0;		// Sto_Rec=1 for wavefiled storage
						// Sto_Rec=0 for wavefiled reconstruction

	int Cut_Sub=0;		// Cut_Sub=1 for Cut direct wave from records
						// Cut_Sub=0 for Substract direct wave by modeling


	int beta1, beta2;				// beta1=0 or 1 for non-disperssive or disperssive; beta2=0 or 1 for non-attenuating or attenuating
	int beta1_c, beta2_c;			// beta1_c=beta1_c for correcting phase distortion; beta2_c=-1*beta2 for compensating amplitude loss

	if(RTMtype==0)			// RTMtype=0 for acoustic RTM
	{
		beta1=0;			// Non-dispersion
		beta2=0;			// Non-attenuation
		beta1_c=0;			// Non-dispersion correction
		beta2_c=0;			// Non-attenuation compensation
		printf("Acoustic RTM is selected!\n");
	}

	else if(RTMtype==1)		// RTMtype=1 for viscoacoustic RTM without compensation
	{
		beta1=1;			// dispersion
		beta2=1;			// attenuation
		beta1_c=0;			// Non-dispersion correction
		beta2_c=0;			// Non-attenuation compensation
		printf("Viscocoustic RTM without compensation is selected!\n");
	}
	else if(RTMtype==2)		// RTMtype=2 for QRTM without stabilization
	{
		beta1=1;			// dispersion
		beta2=1;			// attenuation
		beta1_c=1;			// dispersion correction
		beta2_c=-1;			// attenuation compensation
		printf("Q-RTM without stabilization\n");
	}
	else if(RTMtype==3)		// RTMtype=3 for QRTM using explicit stabilization
	{
		beta1=1;			// dispersion
		beta2=1;			// attenuation
		beta1_c=1;			// dispersion correction
		beta2_c=-1;			// attenuation compensation
		printf("Q-RTM using explicit stabilization alpha=2\n");
	}
	else if(RTMtype==4)		// RTMtype=3 for QRTM using explicit stabilization
	{
		beta1=1;			// dispersion
		beta2=1;			// attenuation
		beta1_c=1;			// dispersion correction
		beta2_c=-1;			// attenuation compensation
		printf("Q-RTM using explicit stabilization alpha=4\n");
	}
	else if(RTMtype==5)		// RTMtype=3 for QRTM using explicit stabilization
	{
		beta1=1;			// dispersion
		beta2=1;			// attenuation
		beta1_c=1;			// dispersion correction
		beta2_c=-1;			// attenuation compensation
		printf("Q-RTM using explicit stabilization alpha=8\n");
	}		
	else
	{
		printf("Please define RTM type!\n");
	}	
											
	int Ckptstype=0;	// Ckptstype=0 for Ave-Distribution
						// Ckptstype=1 for Log-Distribution
						// Ckptstype=2 for Hyb-Distribution 

	int vp_type=0;		// vp_type=0 for homogeneous model
						// vp_type=1 for ture model
						// vp_type=2 for initial model			

	//=========================================================
	//  Parameters of Cartesian coordinate...
	//  =======================================================

	int it,ix,iz,ip,ipp;
	int nz=234;				// vertical samples
	int nx=663;				// horizontal samples
	int L=20;				// Layers of MTF absorbing boundary condition
	int ntz=nz+2*L;			// total vertical samples 
	int ntx=nx+2*L;			// total horizontal samples 
	int nxh=ntx/2;			// half of vertical samples
	int nzh=ntz/2;			// half of horizontal samples
	int ntp=ntz*ntx;		// total samples of the model
	int np=nx*nz;			// total samples of the simulating domain
	float dz=10.0;			// vertical interval
	float dx=10.0;			// horizontal interval
	FILE *fp;				// file pointer
	char filename[40];		// filename char

	/*=========================================================
	  Parameters of ricker wave...
	  ========================================================*/

	float f0=20.0;			// dominate frequency of the ricker
	float t0=1/f0;			// ricker delay time
	float Omega0=20*PI*f0;	// reference frequency
	int nt=2001;			// time samples
	float dt=0.001;			// time interval
	float *ricker;			// ricker wavelet
	ricker=(float *) malloc(nt*sizeof(float));

	if(myid==0)
	{
		ricker_wave(ricker,nt,f0,t0,dt,3); // where parameter 3 means ricker_derivative
		printf("Ricker wave is done!\n");
	}
	MPI_Bcast(ricker,nt,MPI_FLOAT,0,comm);	// broadcast ricker to each node


	//=========================================================
	//  Parameters of Sources and Receivers...
	//  =======================================================

  	int is, irx,irz;
	int nsid,modsr,prcs;
	int iss,eachsid,offsets;

	int ds=16;				// shots interval
	int ns0=18;				// position of the first shot is L+ns0
	int ns=(ntx-2*L-2*ns0)/ds+1;	// total shot number is 4 
	int rnmax=0;					// maxium receiver numbers
	int nrx_obs=100; 				// 100 receivers every side of shot for backward propagation

	struct Source ss[ns];			// struct pointer for source variables
	for(is=0;is<ns;is++)
	{
		ss[is].s_ix=is*ds+L+ns0;	// receiver horizontal position
		ss[is].s_iz=L+5;			// shot vertical position	
		ss[is].r_n=nx; 				// one shot to all receivers			
		ss[is].r_ix=(int*)malloc(sizeof(int)*ss[is].r_n);
		ss[is].r_iz=(int*)malloc(sizeof(int)*ss[is].r_n);
		for(ip=0;ip<ss[is].r_n;ip++)
		{		
			ss[is].r_ix[ip]=L+ip;	// shot horzitonal position
			ss[is].r_iz[ip]=L+5;	// receiver vertical position
		}	
		if(rnmax<ss[is].r_n)
			rnmax=ss[is].r_n;		// maxium receiver numbers
	}
	if(myid==0)
	{
		printf("The total shot number is %d\n",ns);
		printf("The maximum trace number for source is %d\n",rnmax);
	}


	//=========================================================
	//  Parameters of checkpoints...
	//  =======================================================

	int check_steps;			// checkpoints interval
	int N_cp;					// total number of checkpoints
	int *t_cp;					// time point for checkpoints

	if(Ckptstype==0)	// for Ave-Distribution
	{
		check_steps=200;							// every 200 timestep has a checkpoint
		N_cp=2*(int)(nt/check_steps);				// total number of checkpoints (in pair)
		t_cp= (int *) malloc(N_cp*sizeof(int));		// checkpoints pointer
		for(int icp=0;icp<N_cp;icp++)
		{
			if(icp%2==0)
				t_cp[icp]=check_steps*(icp/2+1);	// checkpoints at even timestep
			else
				t_cp[icp]=t_cp[icp-1]+1;			// checkpoints at odd timestep
			if(myid==0)
			{
				printf("checkpoints time is %d\n", t_cp[icp]);
			}
		}			
	}		

	else if(Ckptstype==1)	// Log-Distribution
	{
		check_steps=20;											// every 20*2^n timestep has a checkpoint
		N_cp=2*(int)((log(nt/check_steps)/log(2)+1));

		if(pow(2,N_cp/2-1)+pow(2,N_cp/2-2) < nt/check_steps)
			N_cp+=2;

		t_cp= (int *) malloc(N_cp*sizeof(int));

		for(int icp=0;icp<N_cp;icp++)
		{
			if(icp%2==0)
			{
				t_cp[icp]=check_steps*(int)(pow(2,icp/2));		// checkpoints at even timestep
				if(t_cp[icp]> nt)
					t_cp[icp]-=t_cp[icp-2]/2;
			}
			else
				t_cp[icp]=t_cp[icp-1]+1;						// checkpoints at odd timestep

			if(myid==0)
			{
				printf("checkpoints time is %d\n", t_cp[icp]);
			}		
		}		
	}

	else if(Ckptstype==2)	// Hyb-Distribution
	{
		check_steps=200;
		float min_steps=20;
		int min_N_cp;
		min_N_cp=2*(int)((log(check_steps/min_steps)/log(2)+1));
		if(pow(2,min_N_cp/2-1)+pow(2,min_N_cp/2-2) < check_steps/min_steps)
			min_N_cp+=2;
		N_cp= min_N_cp+2*(int)(nt/check_steps);
		t_cp= (int *) malloc(N_cp*sizeof(int*));

		for(int icp=0;icp<min_N_cp;icp++)
		{
			if(icp%2==0)
			{
				t_cp[icp]=min_steps*(int)(pow(2,icp/2));
				if(t_cp[icp]> check_steps)
					t_cp[icp]-=t_cp[icp-2]/2;
			}
			else
				t_cp[icp]=t_cp[icp-1]+1;

			if(myid==0)
			{
				printf("checkpoints time is %d\n", t_cp[icp]);
			}	
		}
		for(int icp=min_N_cp;icp<N_cp;icp++)
		{
			if(icp%2==0)
				t_cp[icp]=check_steps*((icp-min_N_cp)/2+1);
			else
				t_cp[icp]=t_cp[icp-1]+1;

			if(myid==0)
			{
				printf("checkpoints time is %d\n", t_cp[icp]);
			}	
		}				
	}

	if(myid==0)
	{
		printf("%d Checkpoints are done!\n", N_cp);
	}


	//=========================================================
	//  Parameters of GPU...(we assume each node has the same number of GPUs)
	//  =======================================================

	int i,GPU_N;						// GPU_N stands for the total number of GPUs per node
	getdevice(&GPU_N);					// Obtain the number of GPUs per node
	printf("The available Device number is %d on %s\n",GPU_N,processor_name);
	struct MultiGPU plan[GPU_N];		// struct pointer for MultiGPU variables

	nsid=ns/(GPU_N*numprocs);			// shots number per GPU card
	modsr=ns%(GPU_N*numprocs);			// residual shots after average shot distribution
	prcs=modsr/GPU_N;					// which GPUs at each node will have one more shot
	if(myid<prcs)						
	{
		eachsid=nsid+1;					// if thread ID less than prcs, the corresponding GUPs have one more shot
		offsets=myid*(nsid+1)*GPU_N;	// the offset of the shots
	}
	else
	{
		eachsid=nsid;												// the rest GUPs have nsid shots
		offsets=prcs*(nsid+1)*GPU_N+(myid-prcs)*nsid*GPU_N;			// the offset of the shots (including previous shots)
	}
	

	//=========================================================
	//  Parameters of model...
	//  =======================================================

	float *vp, *Qp;						// velocity and Q models
	float *Gamma, averGamma;			// Q-related parameter Gamma and its average
	float vp_max,Qp_max;				
	float vp_min,Qp_min;
	float avervp;

	vp = (float*)malloc(sizeof(float)*ntp);
	Qp = (float*)malloc(sizeof(float)*ntp);
	Gamma = (float*)malloc(sizeof(float)*ntp);	

	if(myid==0)
	{
		get_acc_model(vp,Qp,ntp,ntx,ntz,L);
		fp=fopen("./output/output_marmousi/acc_vp.dat","wb");
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			for(iz=L;iz<=ntz-L-1;iz++)
			{
				fwrite(&vp[iz*ntx+ix],sizeof(float),1,fp);

			}
		}
		fclose(fp);
		fp=fopen("./output/output_marmousi/acc_Qp.dat","wb");
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			for(iz=L;iz<=ntz-L-1;iz++)
			{
				fwrite(&Qp[iz*ntx+ix],sizeof(float),1,fp);

			}
		}
		fclose(fp);

		vp_max=0.0;
		Qp_max=0.0;
		vp_min=5000.0;
		Qp_min=5000.0;

		for(ip=0;ip<ntp;ip++)
		{     
			if(vp[ip]>=vp_max)
			{
				vp_max=vp[ip];
			}
			if(Qp[ip]>=Qp_max)
			{
				Qp_max=fabs(Qp[ip]);
			}
			if(vp[ip]<=vp_min)
			{
				vp_min=vp[ip];
			}
			if(Qp[ip]<=Qp_min)
			{
				Qp_min=fabs(Qp[ip]);
			}
		}

		printf("vp_max = %f\n",vp_max); 
		printf("Qp_max = %f\n",Qp_max);
		printf("vp_min = %f\n",vp_min); 
		printf("Qp_min = %f\n",Qp_min);

		averGamma=0.0;
		avervp=0.0;

		for(iz=0;iz<=ntz-1;iz++)  
			for(ix=0;ix<=ntx-1;ix++)
			{
				Qp[iz*ntx+ix] /= 1;
				Gamma[iz*ntx+ix]=atan(1/Qp[iz*ntx+ix])/PI;
				averGamma+=Gamma[iz*ntx+ix]/(ntx*ntz);
				avervp+=vp[iz*ntx+ix]/(ntx*ntz);		
			}

		printf("The true model is done!\n"); 
	}

	MPI_Bcast(vp, ntp, MPI_FLOAT, 0, comm);	
	MPI_Bcast(Qp, ntp, MPI_FLOAT, 0, comm);
	MPI_Bcast(&vp_max, 1, MPI_FLOAT, 0, comm);
	MPI_Bcast(&Qp_max, 1, MPI_FLOAT, 0, comm);
	MPI_Bcast(&vp_min, 1, MPI_FLOAT, 0, comm);
	MPI_Bcast(&Qp_min, 1, MPI_FLOAT, 0, comm);	
	MPI_Bcast(Gamma, ntp, MPI_FLOAT, 0, comm);
	MPI_Bcast(&averGamma, 1, MPI_FLOAT, 0, comm);
	MPI_Bcast(&avervp, 1, MPI_FLOAT, 0, comm);

	float *Inner_image_cor, *Inner_image_nor, *Final_image_cor, *Final_image_nor;

	Inner_image_cor=(float*)malloc(sizeof(float)*np);
	Inner_image_nor=(float*)malloc(sizeof(float)*np);
	Final_image_cor=(float*)malloc(sizeof(float)*np);
	Final_image_nor=(float*)malloc(sizeof(float)*np);
	memset(Inner_image_cor,0,np*sizeof(float));
	memset(Inner_image_nor,0,np*sizeof(float));
	memset(Final_image_cor,0,np*sizeof(float));
	memset(Final_image_nor,0,np*sizeof(float));


	// The following two functions are responsible for alloc and initialization struct varibale plan for GPU_N
	cuda_Device_malloc(ntx, ntz, ntp, nx, nz, nt, dx, dz, L, rnmax, N_cp, plan, GPU_N);
	cuda_Host_initialization(ntx, ntz, ntp, nx, nz, nt, dx, dz, L, rnmax, N_cp, plan, GPU_N);





	//=========================================================
	//  Stabilization parameters
	//  =======================================================

	// parameters for explicit satbiliztion scheme
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


	MPI_Barrier(comm);	// MPI barrier to ensure that all variables have been well-defined


	//=======================================================
	//  Calculate the Observed seismograms...
	//  (this process is designed for synthetic example)
	//========================================================

	// This is tradeoff to debug, there is bug in our code when cufft is called first time, it may return random result.
	// Even we have tried to debug over months, this stuble still exists, therefore we add an additional forward simulation.
	// This tradeoff won't affect the proformance of the proposed method. If you could find this bug, please email me.

	is=0;
	vp_type = 1; 
	Save_Not=0;	
	cuda_visco_PSM_2d_forward
	(
		beta1, beta2,
		nt, ntx, ntz, ntp, nx, nz, L, dx, dz, dt,
		vp, Gamma, avervp, averGamma, f0, Omega0, ricker,
		myid, is, ss, plan, GPU_N, rnmax, nrx_obs, N_cp, t_cp, alphaorder, sigmafactor,
		Sto_Rec, vp_type, Save_Not
	);


	for(iss=0; iss<eachsid; iss++)	// each GPU card compute eachside shots 
	{		
		is=offsets+iss*GPU_N;		// current shot index
		vp_type = 1; 				// ture model
		Save_Not=0;					// Don't save any wavefield snapshots

		// viscoacoustic wave equation modeling using ture model to obtain observed seismogram
		cuda_visco_PSM_2d_forward
		(
			beta1, beta2,
			nt, ntx, ntz, ntp, nx, nz, L, dx, dz, dt,
			vp, Gamma, avervp, averGamma, f0, Omega0, ricker,
			myid, is, ss, plan, GPU_N, rnmax, nrx_obs, N_cp, t_cp, alphaorder, sigmafactor,
			Sto_Rec, vp_type, Save_Not
		);

		// write observed seismogram to disk
		for(i=0;i<GPU_N;i++)
		{
			sprintf(filename,"./output/output_marmousi/%dsource_seismogram_obs.dat",is+i+1);
			fp=fopen(filename,"wb");
			for(ix=0;ix<ss[is+i].r_n;ix++)
			{
				for(it=0;it<nt;it++)
				{
					fwrite(&plan[i].seismogram_obs[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
				}
			}
			fclose(fp);
		}

		// remove direct wave from observed seismogram by mutting and obtain the residual seismogram
		if(Cut_Sub==1)
		{
			for(i=0;i<GPU_N;i++)
			{
				cut_dir(plan[i].seismogram_obs, plan[i].seismogram_rms, rnmax, nt, is, dx, dz, dt, 
					ss[is+i].r_iz, ss[is+i].s_ix, ss[is+i].s_iz, t0, vp);
				sprintf(filename,"./input/records_marmousi/%dsource_seismogram_rms_type%d.dat",is+i+1,RTMtype);
				fp=fopen(filename,"wb");
				for(ix=0;ix<ss[is+i].r_n;ix++)
				{
					for(it=0;it<nt;it++)
					{
						fwrite(&plan[i].seismogram_rms[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);
			}				
		}
	}

	if(myid==0)
	{
		printf("seismogram_obs is obtained!\n");
	}
	if(Cut_Sub==1&&myid==0)
	{
		printf("seismogram_rms is obtained!\n");
	}


	//=======================================================
	//  Calculate the direct and rms seismograms...
	//  (this process is designed for removing direct wave from 
	//  observed seismogram by subtracting direct wave, which is
	//  simulated by homogeneous velocity and Q model)
	//========================================================

	if(Cut_Sub==0)
	{
		if(myid==0)
		{
			get_homo_model(vp,ntp,ntx,ntz,L);
			get_homo_model(Gamma,ntp,ntx,ntz,L);
			fp=fopen("./output/output_marmousi/homo_vp.dat","wb");
			for(ix=L;ix<=ntx-L-1;ix++)
			{
				for(iz=L;iz<=ntz-L-1;iz++)
				{
					fwrite(&vp[iz*ntx+ix],sizeof(float),1,fp);
				}
			}
			fclose(fp);
			printf("The homogeneous model is done!\n"); 			
		}
		MPI_Bcast(vp, ntp, MPI_FLOAT, 0, comm);
		MPI_Bcast(Gamma, ntp, MPI_FLOAT, 0, comm);

		for(int iss=0; iss<eachsid; iss++)  
		{
			is=offsets+iss*GPU_N;
			vp_type = 0; // homogeneous model
			Save_Not=0;

			//// viscoacoustic wave equation modeling using homogeneous model to obtain direct seismogram
			cuda_visco_PSM_2d_forward
			(
				beta1, beta2,
				nt, ntx, ntz, ntp, nx, nz, L, dx, dz, dt,
				vp, Gamma, avervp, averGamma, f0, Omega0, ricker,
				myid, is, ss, plan, GPU_N, rnmax, nrx_obs, N_cp, t_cp, alphaorder, sigmafactor,
				Sto_Rec, vp_type, Save_Not
			);

			// write direct seismogram to disk and calculate residual seismogram 
			for(i=0;i<GPU_N;i++)
			{
				sprintf(filename,"./output/output_marmousi/%dsource_seismogram_dir.dat",is+i+1);
				fp=fopen(filename,"wb");
				for(ix=0;ix<ss[is+i].r_n;ix++)
				{
					for(it=0;it<nt;it++)
					{
						fwrite(&plan[i].seismogram_dir[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);

				sprintf(filename,"./output/output_marmousi/%dsource_seismogram_obs.dat",is+i+1);
				fp=fopen(filename,"rb");
				for(ix=0;ix<ss[is+i].r_n;ix++)
				{
					for(it=0;it<nt;it++)
					{
						fread(&plan[i].seismogram_obs[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);

				for(ip=0;ip<ss[is+i].r_n*nt;ip++)
				{
					plan[i].seismogram_rms[ip]=plan[i].seismogram_obs[ip]-plan[i].seismogram_dir[ip];
				}

				sprintf(filename,"./input/records_marmousi/%dsource_seismogram_rms_type%d.dat",is+i+1,RTMtype);
				fp=fopen(filename,"wb");
				for(ix=0;ix<ss[is+i].r_n;ix++)
				{
					for(it=0;it<nt;it++)
					{
						fwrite(&plan[i].seismogram_rms[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);
			}
		}
		if(myid==0)
		{
			printf("seismogram_dir is obtained!\n");
			printf("seismogram_rms is obtained!\n");
		}
	}



	//=======================================================
	//  Construct the forward wavefields and Back-propagate
	//  the RMS seismograms, Meanwhile the images are computed... 
	//========================================================

	if(myid==0)
	{
		printf("====================\n");
		printf("    RTM BEGIN\n");
		printf("====================\n");

		start=clock();				// start clock
	}

	// obtain the initial model for forward and backward propagation
	if(myid==0)
	{
		// obtain the ture velocity and Q model
		get_acc_model(vp,Qp,ntp,ntx,ntz,L);
		for(iz=0;iz<=ntz-1;iz++)  
			for(ix=0;ix<=ntx-1;ix++)
			{
				Qp[iz*ntx+ix] /= 1;
				Gamma[iz*ntx+ix]=atan(1/Qp[iz*ntx+ix])/PI;	
			}
		printf("The true model is done!\n"); 

		// obtain the initial model
		get_ini_model(vp,ntp,ntx,ntz,L);
		fp=fopen("./output/output_marmousi/ini_vp.dat","wb");
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			for(iz=L;iz<=ntz-L-1;iz++)
			{
				fwrite(&vp[iz*ntx+ix],sizeof(float),1,fp);
			}
		}
		fclose(fp);
		printf("The initial model is done!\n"); 
	}
	MPI_Bcast(vp, ntp, MPI_FLOAT, 0, comm);
	MPI_Bcast(Gamma, ntp, MPI_FLOAT, 0, comm);


	// debug part
	is=0;
	vp_type = 2;				
	Save_Not=0;
	cuda_visco_PSM_2d_forward
	(
		beta1_c, beta2_c,
		nt, ntx, ntz, ntp, nx, nz, L, dx, dz, dt,
		vp, Gamma, avervp, averGamma, f0, Omega0, ricker,
		myid, is, ss, plan, GPU_N, rnmax, nrx_obs, N_cp, t_cp, alphaorder, sigmafactor,
		Sto_Rec, vp_type, Save_Not
	);


	for(int iss=0; iss<eachsid; iss++)  
	{
		is=offsets+iss*GPU_N;
		vp_type = 2;				// initial model
		Save_Not=0;

		// compensated source wavefiled propagation using initial velocity
		cuda_visco_PSM_2d_forward
		(
			beta1_c, beta2_c,
			nt, ntx, ntz, ntp, nx, nz, L, dx, dz, dt,
			vp, Gamma, avervp, averGamma, f0, Omega0, ricker,
			myid, is, ss, plan, GPU_N, rnmax, nrx_obs, N_cp, t_cp, alphaorder, sigmafactor,
			Sto_Rec, vp_type, Save_Not
		);

		// write the synthetic seismogram to disk
		for(i=0;i<GPU_N;i++)
		{
			sprintf(filename,"./output/output_marmousi/%dsource_seismogram_syn.dat",is+i+1);
			fp=fopen(filename,"wb");
			for(ix=0;ix<ss[is+i].r_n;ix++)
			{
				for(it=0;it<nt;it++)
				{
					fwrite(&plan[i].seismogram_syn[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
				}
			}
			fclose(fp);
		}  

		// read the residual seismogram from disk
		for(i=0;i<GPU_N;i++)
		{
			sprintf(filename,"./input/records_marmousi/%dsource_seismogram_rms_type%d.dat",is+i+1,RTMtype);
			fp=fopen(filename,"rb");
			for(ix=0;ix<ss[is+i].r_n;ix++)
			{
				for(it=0;it<nt;it++)
				{
					fread(&plan[i].seismogram_rms[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
				}
			}
			fclose(fp);
		}

		// compensated receiver wavefiled propagation and imaging using initial velocity
		cuda_visco_PSM_2d_backward
		(
			Geometry, beta1_c, beta2_c,
			nt, ntx, ntz, ntp, nx, nz, L, dx, dz, dt,
			vp, Gamma, avervp, averGamma, f0, Omega0, ricker,
			myid, is, ss, plan, GPU_N, rnmax, nrx_obs, N_cp, t_cp, alphaorder, sigmafactor,
			Sto_Rec, Save_Not
		);

		// output images
		for(i=0;i<GPU_N;i++)
		{
			
			// calculate normalized image
			float image_sources_max=0.0;
			for(ip=0;ip<ntp;ip++)
			{
				if(image_sources_max<fabs(plan[i].image_sources[ip]))
				{
					image_sources_max=fabs(plan[i].image_sources[ip]);
				}
			}
			for(ip=0;ip<ntp;ip++)		
			{
				plan[i].image_nor[ip]=plan[i].image_cor[ip]
						/(plan[i].image_sources[ip]+1.0e-5*image_sources_max);
			}

			// Laplace filtering
			Laplace_FD_filtering(plan[i].image_cor,ntx,ntz,dx,dz);
			Laplace_FD_filtering(plan[i].image_nor,ntx,ntz,dx,dz);

			// write all images to disk
			sprintf(filename,"./output/output_marmousi/image_sources%d.dat",is+i+1);
			fp=fopen(filename,"wb");
			fwrite(&plan[i].image_sources[0],sizeof(float),ntp,fp);
			fclose(fp);
			sprintf(filename,"./output/output_marmousi/image_receivers%d.dat",is+i+1);
			fp=fopen(filename,"wb");
			fwrite(&plan[i].image_receivers[0],sizeof(float),ntp,fp);
			fclose(fp);
			sprintf(filename,"./output/output_marmousi/image_cor%d.dat",is+i+1);
			fp=fopen(filename,"wb");
			fwrite(&plan[i].image_cor[0],sizeof(float),ntp,fp);
			fclose(fp);
			sprintf(filename,"./output/output_marmousi/image_nor%d.dat",is+i+1);
			fp=fopen(filename,"wb");
			fwrite(&plan[i].image_nor[0],sizeof(float),ntp,fp);
			fclose(fp);


			sprintf(filename,"./output/output_marmousi/image_cor%d.dat",is+i+1);
			fp=fopen(filename,"rb");
			fread(&plan[i].image_cor[0],sizeof(float),ntp,fp);
			fclose(fp);

			sprintf(filename,"./output/output_marmousi/image_nor%d.dat",is+i+1);
			fp=fopen(filename,"rb");
			fread(&plan[i].image_nor[0],sizeof(float),ntp,fp);
			fclose(fp);

			// calculate inner images at simulation domain (exclude absorbing boundary) 
			for(iz=L;iz<=ntz-L-1;iz++)
			{
				for(ix=L;ix<=ntx-L-1;ix++)
				{
					ip=iz*ntx+ix;
					ipp=(ix-L)*nz+iz-L;
					Inner_image_cor[ipp]+=plan[i].image_cor[ip]*vp[ip]*vp[ip];
					Inner_image_nor[ipp]+=plan[i].image_nor[ip]*vp[ip]*vp[ip];
				}
			}
		}
	}//end is (shotnumbers)

	// MPI barrier to ensure all shots have been calculted before stacking
	MPI_Barrier(comm);

	// stacking all shot's images
	MPI_Allreduce(Inner_image_cor,Final_image_cor,np,MPI_FLOAT,MPI_SUM,comm);
	MPI_Allreduce(Inner_image_nor,Final_image_nor,np,MPI_FLOAT,MPI_SUM,comm);

	//==========================================================
	//  Output the final images,...
	//===========================================================
	if(myid==0)
	{
		sprintf(filename,"./output/output_marmousi/Final_image_cor_type%d.dat",RTMtype);
		fp=fopen(filename,"wb");
		fwrite(&Final_image_cor[0],sizeof(float),np,fp);
		fclose(fp);
		sprintf(filename,"./output/output_marmousi/Final_image_nor_type%d.dat",RTMtype);
		fp=fopen(filename,"wb");
		fwrite(&Final_image_nor[0],sizeof(float),np,fp);
		fclose(fp);
	}

	MPI_Barrier(comm);

	if(myid==0)
	{
		printf("====================\n");
		printf("      THE END\n");
		printf("====================\n");

		end=clock();
		printf("The cost of the run time for Marmousi expQRTM is %f seconds\n",
				(double)(end-start)/CLOCKS_PER_SEC);
	}

	//==========================================================
	//  Free the variables...
	//==========================================================

	// the following function is responsible for free struct variable plan for GPU_N
	cuda_Device_free(ntx, ntz, ntp, nx, nz, nt, dx, dz, L, rnmax, N_cp, plan, GPU_N);

	for(is=0;is<ns;is++)
	{
		free(ss[is].r_ix);
		free(ss[is].r_iz);
	} 
	free(ricker);
	free(vp); free(Qp); free(Gamma);
	free(t_cp);
	free(Inner_image_cor);
	free(Inner_image_nor);
	free(Final_image_cor);
	free(Final_image_nor);

	// MPI end
	MPI_Barrier(comm);
	MPI_Finalize();
	
	return 0;
}






//==========================================================
//  This subroutine is used for calculating the ricker wave
//  ========================================================
void ricker_wave
(
	float *ricker, int nt, float f0, float t0, float dt, int flag
)
{
	float pi=3.1415927;
	int   it;
	float temp,max=0.0;
	FILE *fp;
	if(flag==1)
	{
		for(it=0;it<nt;it++)
		{
			temp=pi*f0*(it*dt-t0);
			temp=temp*temp;
			ricker[it]=(1.0-2.0*temp)*exp(-temp);
		}
		fp=fopen("./output/output_marmousi/ricker.dat","wb");    
		for(it=0;it<nt;it++)
		{
			fwrite(&ricker[it],sizeof(float),1,fp);
		}    
		fclose(fp);
	}
	if(flag==2)
	{
		for(it=0;it<nt;it++)
		{
			temp=pi*f0*(it*dt-t0);
			temp=temp*temp;         
			ricker[it]=(it*dt-t0)*exp(-temp);

			if(max<fabs(ricker[it]))
			{
				max=fabs(ricker[it]);
			}
		}
		for(it=0;it<nt;it++)
		{
			ricker[it]=ricker[it]/max;
		}
		fp=fopen("./output/output_marmousi/ricker_integration.dat","wb");    
		for(it=0;it<nt;it++)
		{
			fwrite(&ricker[it],sizeof(float),1,fp);
		}    
		fclose(fp);
	}
	if(flag==3)
	{	
		for(it=0;it<nt;it++)
		{
			temp=pi*f0*(it*dt-t0);
			ricker[it]=(4*powf(pi*f0,4)*powf((it*dt-t0),3)-6*powf(pi*f0,2)*(it*dt-t0))*exp(-powf(temp,2));  

			if(max<fabs(ricker[it]))
			{
				max=fabs(ricker[it]);
			}
		}
		for(it=0;it<nt;it++)
		{
			ricker[it]=ricker[it]/max;
		}
		fp=fopen("./output/output_marmousi/ricker_derivative.dat","wb");    
		for(it=0;it<nt;it++)
		{
			fwrite(&ricker[it],sizeof(float),1,fp);
		}    
		fclose(fp);
	}
	return;
}


/*==========================================================
  This subroutine is used for initializing the true model...
  ===========================================================*/
void get_acc_model
(
	float *vp, float *Qp, int ntp, int ntx, int ntz, int L
)
{
	int ip,ipp,iz,ix; 
	FILE *fp;
	fp=fopen("./input/para_marmousi/acc_vp.dat","rb");
	for(ix=L;ix<ntx-L;ix++)
	{
		for(iz=L;iz<ntz-L;iz++)
		{
			ip=iz*ntx+ix;
			fread(&vp[ip],sizeof(float),1,fp);           
		}
	}
	fclose(fp);
	for(iz=0;iz<=L-1;iz++)
	{
		for(ix=0;ix<=L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=L*ntx+L;

			vp[ip]=vp[ipp];
		}
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=L*ntx+ix;

			vp[ip]=vp[ipp];
		}
		for(ix=ntx-L;ix<ntx;ix++)
		{
			ip=iz*ntx+ix;
			ipp=L*ntx+ntx-L-1;

			vp[ip]=vp[ipp];
		}
	}
	for(iz=L;iz<=ntz-L-1;iz++)
	{
		for(ix=0;ix<=L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=iz*ntx+L;

			vp[ip]=vp[ipp];
		}
		for(ix=ntx-L;ix<ntx;ix++)
		{
			ip=iz*ntx+ix;
			ipp=iz*ntx+ntx-L-1;

			vp[ip]=vp[ipp];
		}
	}

	for(iz=ntz-L;iz<ntz;iz++)
	{
		for(ix=0;ix<=L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=(ntz-L-1)*ntx+L;

			vp[ip]=vp[ipp];
		}
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=(ntz-L-1)*ntx+ix;

			vp[ip]=vp[ipp];
		}
		for(ix=ntx-L;ix<ntx;ix++)
		{
			ip=iz*ntx+ix;
			ipp=(ntz-L-1)*ntx+ntx-L-1;
			vp[ip]=vp[ipp];
		}
	}
	fp=fopen("./input/para_marmousi/acc_Qp.dat","rb");
	for(ix=L;ix<ntx-L;ix++)
	{
		for(iz=L;iz<ntz-L;iz++)
		{
			ip=iz*ntx+ix;
			fread(&Qp[ip],sizeof(float),1,fp);
			Qp[ip]=Qp[ip];
		}
	}
	fclose(fp);
	for(iz=0;iz<=L-1;iz++)
	{
		for(ix=0;ix<=L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=L*ntx+L;
			Qp[ip]=Qp[ipp];
		}
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=L*ntx+ix;
			Qp[ip]=Qp[ipp];
		}
		for(ix=ntx-L;ix<ntx;ix++)
		{
			ip=iz*ntx+ix;
			ipp=L*ntx+ntx-L-1;
			Qp[ip]=Qp[ipp];
		}
	}
	for(iz=L;iz<=ntz-L-1;iz++)
	{
		for(ix=0;ix<=L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=iz*ntx+L;
			Qp[ip]=Qp[ipp];
		}
		for(ix=ntx-L;ix<ntx;ix++)
		{
			ip=iz*ntx+ix;
			ipp=iz*ntx+ntx-L-1;
			Qp[ip]=Qp[ipp];
		}
	}

	for(iz=ntz-L;iz<ntz;iz++)
	{
		for(ix=0;ix<=L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=(ntz-L-1)*ntx+L;
			Qp[ip]=Qp[ipp];
		}
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=(ntz-L-1)*ntx+ix;
			Qp[ip]=Qp[ipp];
		}
		for(ix=ntx-L;ix<ntx;ix++)
		{
			ip=iz*ntx+ix;
			ipp=(ntz-L-1)*ntx+ntx-L-1;
			Qp[ip]=Qp[ipp];
		}
	}
	return;
}


// ==========================================================
//  This subroutine is used for initializing the homogeneous model...
//  =========================================================
void get_homo_model
(
	float *vp, int ntp, int ntx, int ntz, int L
)
{
	int ip,ipp,iz,ix;  
	FILE *fp;
	for(ix=0;ix<ntx;ix++)
	{
		for(iz=0;iz<ntz;iz++)
		{
			ip=iz*ntx+ix;
			if(iz>L+1)
			{
				ipp=(L+1)*ntx+ix;
				vp[ip]=vp[ipp];
			}
		}
	}
	return;
}


//==========================================================
//  This subroutine is used for initializing the initial model...
//  ========================================================
void get_ini_model
(
	float *vp, int ntp, int ntx, int ntz, int span
)
{
	int ix, ixw, ixx;
	int iz, izw, izz;

	float *s_a;
	s_a=(float*)malloc(sizeof(float)*ntp);

	float *sx_a;
	sx_a=(float*)malloc(sizeof(float)*ntp);

	for(iz=0;iz<ntz;iz++)
	{
		for(ix=0; ix<ntx; ix++)
		{
			s_a[iz*ntx+ix]=0.0;
			for(ixw=ix-span;ixw<=ix+span;ixw++)
			{
				if(ixw<0)
					ixx=0;
				else if(ixw>ntx-1)
					ixx=ntx-1;
				else
					ixx=ixw;
				s_a[iz*ntx+ix]+=vp[iz*ntx+ixx]/(2*span+1);
			}		
		}	
	}

	for(iz=0;iz<ntz;iz++)
	{
		for(ix=0; ix<ntx; ix++)
		{		
			sx_a[iz*ntx+ix]=s_a[iz*ntx+ix];		
		}	
	}

	for(iz=0;iz<ntz;iz++)
	{
		for(ix=0; ix<ntx; ix++)
		{
			s_a[iz*ntx+ix]=0.0;
			for(izw=iz-span;izw<=iz+span;izw++)
			{
				if(izw<0)
					izz=0;
				else if(izw>ntz-1)
					izz=ntz-1;
				else
					izz=izw;
				s_a[iz*ntx+ix]+=sx_a[izz*ntx+ix]/(2*span+1);
			}		
		}	
	}

	for(iz=0;iz<ntz;iz++)
	{
		for(ix=0; ix<ntx; ix++)
		{
			vp[iz*ntx+ix]=s_a[iz*ntx+ix];		
		}	
	}

	free(s_a);	
	free(sx_a);		
}


//==========================================================
//  This subroutine is used for cutting the direct wave
//  =======================================================

void cut_dir
(
	float *seismogram_obs, float *seismogram_rms, 
	int rnmax, int nt, int is, float dx, float dz, float dt, 
	int *r_iz, int s_ix, int s_iz, float t0, float *vp
)
{
	int it, ix;
	float removeline[rnmax];

	for(it=0;it<nt; it++)
	{
		for(ix=0; ix<rnmax; ix++)
		{
			removeline[ix]=(sqrt(powf((ix-s_ix)*dx,2)+powf((r_iz[0]-s_iz)*dz,2))/vp[1*rnmax+ix]+4.0*t0)/dt;

			if(it<removeline[ix])
				seismogram_rms[it*rnmax+ix]=0.0;
			else
				seismogram_rms[it*rnmax+ix]=seismogram_obs[it*rnmax+ix];
		}	
	}
}

//==========================================================
//  This subroutine is used for Laplace filtering
//  ========================================================
void Laplace_FD_filtering
(
	float *image, int ntx, int ntz, float dx, float dz
)
{ 
	int ix,iz,ip;
	float diff1, diff2;
	float *tmp;
	tmp = (float*)malloc(sizeof(float)*ntx*ntz);
	memset(tmp, 0, ntx*ntz*sizeof(float));
	for(iz=1;iz<ntz-1;iz++)
	{
		for(ix=1;ix<ntx-1;ix++)
		{
			ip=iz*ntx+ix;
			diff1=(image[ip+ntx]-2.0*image[ip]+image[ip-ntx])/(dz*dz);
			diff2=(image[ip+1]-2.0*image[ip]+image[ip-1])/(dx*dx);	
			tmp[ip]=diff1+diff2;          
		}
	}

	for(iz=0;iz<=ntz-1;iz++)
	{
		for(ix=0;ix<=ntx-1;ix++)
		{
			ip=iz*ntx+ix;
			image[ip]=tmp[ip];          
		}
	}	
	free(tmp);
	return;
}
