// This is the header file for QRTM code package. We define two struct variables in this header 
// and put all function statement in here for both extern "C" functions in CUDAQRTM.cu and C 
// functions in CQRTM.cpp.


#include "cufft.h"

extern "C"
void getdevice(int *GPU_N);

extern "C"
struct Source
{
	int s_iz, s_ix, *r_iz, *r_ix, r_n;
};



//=======================================================
//  Defining the struct variable which contains all 
//	variables allocated in device.
//=======================================================

extern "C"
struct MultiGPU
{
	// cufft handle for forward and backward propagation
	cufftHandle PLAN_FORWARD;
	cufftHandle PLAN_BACKWARD;

	// host pagelock memory (variables needs to cudaMemcpyDeviceToHost)
	cufftComplex *u0, *u1, *u2;
	float *seismogram_obs;
	float *seismogram_dir;
	float *seismogram_syn;
	float *seismogram_rms;
	float *image_sources, *image_receivers;
	float *image_cor, *image_nor;

	// device global memory (variables needs to cudaMemcpyHostToDevice)
	int *d_r_ix;
	int *d_r_iz;
	float *d_ricker;
	float *d_vp, *d_Gamma;
	cufftComplex *d_u0, *d_u1, *d_u2;
	cufftComplex *d_u0_inv, *d_u1_inv, *d_u2_inv;
	int *d_t_cp;
	float *d_u_cp;
	float *d_kx, *d_kz;
	cufftComplex *d_uk, *d_uk0;
	cufftComplex *d_uk_inv, *d_uk0_inv;
	cufftComplex *d_Lap_uk, *d_amp_uk, *d_pha_uk, *d_sta_uk;
	cufftComplex *d_Lap, *d_amp_Lap, *d_pha_Lap, *d_sta_Lap;
	float *d_seismogram;
	float *d_seismogram_rms;
	float *d_borders_up,*d_borders_bottom;
	float *d_borders_left,*d_borders_right;
	float *d_u2_final0, *d_u2_final1;
	float *d_image_sources,*d_image_receivers;
	float *d_image_cor, *d_image_nor;
};



//=======================================================
//  Statements for all C function in QRTM_*.cpp 
//=======================================================
void ricker_wave
(
	float *ricker, int nt, float f0, float t0, float dt, int flag
);

void get_acc_model
(
	float *vp, float *Qp, int ntp, int ntx, int ntz, int L
);

void get_homo_model
(
	float *vp, int ntp, int ntx, int ntz, int L
);

void get_ini_model
(
	float *vp, int ntp, int ntx, int ntz, int span
);

void cut_dir
(
	float *seismogram_obs, float *seismogram_rms, 
	int rnmax, int nt, int is, float dx, float dz, float dt, 
	int *r_iz, int s_ix, int s_iz, float t0, float *vp
);

void Laplace_FD_filtering
(
	float *image, int ntx, int ntz, float dx, float dz
);


//=======================================================
//  Statements for all extern "C" function in CUDAQRTM.cu
//=======================================================
extern "C"
void cuda_visco_PSM_2d_forward
(
	int beta1, int beta2,
	int nt, int ntx, int ntz, int ntp, int nx, int nz, int L, float dx, float dz, float dt,
	float *vp, float *Gamma, float avervp, float averGamma, float f0, float Omega0, float *ricker,
	int myid, int is, struct Source ss[], struct MultiGPU plan[], int GPU_N, int rnmax, int nrx_obs, int N_cp, int *t_cp,
	float alphaorder, float sigmafactor,
	int Sto_Rec, int vp_type, int Save_Not
);

extern "C"
void cuda_visco_PSM_2d_backward
(
	int Geometry, int beta1, int beta2,
	int nt, int ntx, int ntz, int ntp, int nx, int nz, int L, float dx, float dz, float dt,
	float *vp, float *Gamma, float avervp, float averGamma, float f0, float Omega0, float *ricker,
	int myid, int is, struct Source ss[], struct MultiGPU plan[], int GPU_N, int rnmax, int nrx_obs, int N_cp, int *t_cp,
	float alphaorder, float sigmafactor,
	int Sto_Rec, int Save_Not
);

extern "C"
void Laplace_filtering
(
	float *image, int ntx, int ntz, float dx, float dz
);

extern "C"
void cuda_Device_malloc
(
	int ntx, int ntz, int ntp, int nx, int nz, int nt, 
	float dx, float dz, int L, int rnmax, int N_cp,
	struct MultiGPU plan[], int GPU_N
);

extern "C"
void cuda_Device_free
(
	int ntx, int ntz, int ntp, int nx, int nz, int nt, 
	float dx, float dz, int L, int rnmax, int N_cp,
	struct MultiGPU plan[], int GPU_N
);

extern "C"
void cuda_Host_initialization
(
	int ntx, int ntz, int ntp, int nx, int nz, int nt, 
	float dx, float dz, int L, int rnmax, int N_cp,
	struct MultiGPU plan[], int GPU_N
);



