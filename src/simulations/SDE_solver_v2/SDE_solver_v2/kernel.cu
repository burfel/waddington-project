#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda.h> 
#include <curand.h> 

#define CUDA_CALL(x) do { if((x)!=cudaSuccess){ \
	printf("Error at %s:%d\n",__FILE__,__LINE__); \
	return EXIT_FAILURE;}} while(0) 

#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
	printf("Error at %s:%d\n",__FILE__,__LINE__); \
	return EXIT_FAILURE;}} while(0)

using namespace std;



__global__ void SDE_simulations(double *x1, double *x2, double *temp1, double *temp2, int n, int time, float *d_Normal, float *d_Uni)
{
	int i = threadIdx.x;
	double a1 = 1.0;
	double a2 = 1.0;
	double b1 = 1.0;
	double b2 = 1.0;
	double k1 = 1.0;
	double k2 = 1.0;
	double H = 4.0;
	double S = 0.5;

	double dt = 10.0 / (double)time;


	if (i<n)
	{
		x1[i] = (double)d_Uni[2*i]*3.0;
		x2[i] = (double)d_Uni[2*i+1]*3.0;
		//x1[i] = 3.0;
		//x2[i] = 3.0;

		for (int j = 0; j < time - 1; j++)
		{
			temp1[i] = x1[i];
			temp2[i] = x2[i];

			x1[i] = temp1[i] + dt*(a1*pow(temp1[i], H) / (pow(temp1[i], H) + pow(S, H)) + b1*pow(S, H) / (pow(temp2[i], H) + pow(S, H)) - k1*temp1[i])+ sqrt(fabs(a1*pow(temp1[i], H) / (pow(temp1[i], H) + pow(S, H)) + b1*pow(S, H) / (pow(temp2[i], H) + pow(S, H)) - k1*temp1[i]))*(double)d_Normal[i*time + 2 * j];
			x2[i] = temp2[i] + dt*(a2*pow(temp2[i], H) / (pow(temp2[i], H) + pow(S, H)) + b2*pow(S, H) / (pow(temp1[i], H) + pow(S, H)) - k2*temp2[i])+ sqrt(fabs(a2*pow(temp2[i], H) / (pow(temp2[i], H) + pow(S, H)) + b2*pow(S, H) / (pow(temp1[i], H) + pow(S, H)) - k2*temp2[i]))*(double)d_Normal[i*time + 2 * j + 1];
		}
	}
}



int main()
{
	int n = 1000000;
	int time = 10;

	//CPU Memory variables
	double *x1, *x2 , *temp1, *temp2;
	x1 = (double *)malloc(n* sizeof(double));
	x2 = (double *)malloc(n* sizeof(double));
	temp1 = (double *)malloc(n * sizeof(double));
	temp2 = (double *)malloc(n * sizeof(double));

	//GPU Memory variables
	double *d_x2, *d_x1 , *d_temp1, *d_temp2;

	cudaMalloc(&d_x1, n * sizeof(double));
	cudaMalloc(&d_x2, n * sizeof(double));
	cudaMalloc(&d_temp1, n * sizeof(double));
	cudaMalloc(&d_temp2, n * sizeof(double));

	curandGenerator_t gen;
	float *d_Normal, *d_Uni;//, *hostData;
	//hostData = (double *)calloc(n, sizeof(double));

	CUDA_CALL(cudaMalloc(&d_Normal, 2 * time*n * sizeof(float)));
	CUDA_CALL(cudaMalloc(&d_Uni, 2*n * sizeof(float)));
	CURAND_CALL(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT));
	CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen, 1234ULL));
	CURAND_CALL(curandGenerateNormal(gen, d_Normal, 2 * n*time, 0.0, 10.0/(float)time));
	CURAND_CALL(curandGenerateUniform(gen, d_Uni, 2 * n));
	

	/*for (int i = 0; i < n; i++)
	{
		x1[i] = 1.0;
		x2[i] = 1.0;
	}*/

	//cudaMemcpy(d_x1, x1, n* sizeof(double), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_x2, x2, n* sizeof(double), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_temp1, temp1, n*sizeof(double), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_temp2, temp2, n*sizeof(double), cudaMemcpyHostToDevice);

	SDE_simulations<<<1,512>>>(d_x1, d_x2, d_temp1, d_temp2, n, time, d_Normal, d_Uni);

	cudaMemcpy(x1, d_x1, n*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(x2, d_x2, n*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(temp1, d_temp1, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(temp2, d_temp2, n * sizeof(double), cudaMemcpyDeviceToHost);

	for (int i = 0; i < 10; i++)
	{
		printf("%f - %f\n", x1[i], x2[i]);
	}


	CURAND_CALL(curandDestroyGenerator(gen));
	CUDA_CALL(cudaFree(d_Normal));
	CUDA_CALL(cudaFree(d_Uni));
	free(x1);
	free(x2);
	free(temp1);
	free(temp2);
	cudaFree(d_x1);
	cudaFree(d_x2);
	cudaFree(d_temp1);
	cudaFree(d_temp2);

	return 0;
}

