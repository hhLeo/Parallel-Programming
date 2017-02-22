//kn.cu
//nvcc -arch=sm_30 kn.cu -o kn.out
/*
时间/s	求距离	求knn	总计
Text 1024	0.003796	0.004046	0.007852
Text 4096	0.231635	0.025109	0.256753
Text 16384	14.722307	0.410728	15.133045
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <fstream>
#include <time.h>
using namespace std;

#define MAXM 16384
#define MAXN 4096
#define MAXK 32
#define MAX_VALUE 5e5	//20^2*1000 = 4e5
#define BLOCK 32 		//max = 32

int m, n, k, *element, *mul, *result;

/* seq */
int getSum(int i, int j) {
	int tmpSum = 0;
	for (int g = 0; g < n; g++) {
		tmpSum += ((element[i*n+g] - element[j*n+g]) * (element[i*n+g] - element[j*n+g]));
	}
	return tmpSum;
}

__global__ void getSum1(int m, int n, int* d_element, int* d_mul) {
	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
		
	int i = BLOCK * by + ty;
	int j = BLOCK * bx + tx;

	int tmpSum = 0;

	__shared__ int SA[BLOCK][BLOCK];
	__shared__ int SB[BLOCK][BLOCK];

	for (int g = 0; g < n; g += BLOCK) {
		//parallel load
		SA[ty][tx] = d_element[i*n+(g+tx)];
		SB[ty][tx] = d_element[j*n+(g+ty)];

		__syncthreads();

		//compute submatrix
		for(int w = 0; w < BLOCK; w++) {
			int tmp = SA[ty][w] - SB[w][tx];
			tmpSum += tmp * tmp;
		}
		__syncthreads();
		//output
		d_mul[i*m+j] = tmpSum;
	}
}

//avoid the duplicated computation
__global__ void getSum2(int m, int n, int* d_element, int* d_mul) {
	int bx = blockIdx.x;
	int by = blockIdx.y;
	if (bx > by)	return;

	int tx = threadIdx.x;
	int ty = threadIdx.y;
		
	int i = BLOCK * by + ty;
	int j = BLOCK * bx + tx;

	int tmpSum = 0;

	__shared__ int SA[BLOCK][BLOCK];
	__shared__ int SB[BLOCK][BLOCK];

	for (int g = 0; g < n; g += BLOCK) {
		//parallel load
		SA[ty][tx] = d_element[i*n+(g+tx)];
		SB[ty][tx] = d_element[j*n+(g+ty)];

		__syncthreads();

		//compute submatrix
		for(int w = 0; w < BLOCK; w++) {
			int tmp = SA[ty][w] - SB[w][tx];
			tmpSum += tmp * tmp;
		}
		__syncthreads();
		//output
		if (i == j)	tmpSum = MAX_VALUE;
		d_mul[j*m+i] = tmpSum;
		if (bx <= by)	d_mul[i*m+j] = tmpSum;
	}
}

//seq
void getknn(int tt) {
	int i = tt;
	int tmpMul[32];

	for (int j = 0; j < k; j++) {
		tmpMul[j] = MAX_VALUE;
	}

	for (int j = 0; j < m; j++) {
		int g;
		for(g = k-1; g >= 0; g--) {
			if(mul[i*m+j] >= tmpMul[g]) {
				break;
			}
		}
		if(g < k - 1) {
			for(int v = k - 2; v > g; v--) {
				result[i*k+v+1] = result[i*k+v];
				tmpMul[v+1] = tmpMul[v];
			}
			result[i*k+g+1] = j;
			tmpMul[g+1] = mul[i*m+j];
		}
	}
}

__global__ void getknn5(int m,int k,int* mul,int* result) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	int tmpMul[32];

	for (int j = 0; j < k; j++) {
		tmpMul[j] = MAX_VALUE;
	}

	for (int j = 0; j < m; j++) {
		int g;
		for(g = k-1; g >= 0; g--) {
			if(mul[i*m+j] >= tmpMul[g]) {
				break;
			}
		}
		if(g < k - 1) {
			for(int v = k - 2; v > g; v--) {
				result[i*k+v+1] = result[i*k+v];
				tmpMul[v+1] = tmpMul[v];
			}
			result[i*k+g+1] = j;
			tmpMul[g+1] = mul[i*m+j];
		}
	}
}

int main(int argc, char const *argv[]) {
	cudaSetDevice(1);
	if (argc != 3)
	{
		printf("usage: ./a.out <input filename> <output filename>\n");
		return 0;
	}

	/* input */
	ifstream fin;
	fin.open(argv[1]);
	ofstream fout;
	fout.open(argv[2]);
	element = (int*)malloc(sizeof(int)*MAXM*MAXN);
	mul = (int*)malloc(sizeof(int)*MAXM*MAXM);
	result = (int*)malloc(sizeof(int)*MAXM*MAXK);
	fin >> m >> n >> k;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			fin >> element[i*n+j];

	if (m < 1024) {
		for (int i = 0; i < m; i++) {
			mul[i*m+i] = MAX_VALUE;
			for (int j = i+1; j < m; j++) {
				mul[i*m+j] = getSum(i, j);
				mul[j*m+i] = mul[i*m+j];
			}
		}
		for (int i = 0; i < m; i++) {
			getknn(i);
		}
		/* output */
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < k; j++) {
				printf("%d ", result[i*k+j]);
				fout << result[i*k+j] << ' ';
			}
			printf("\n");
			fout << endl;
		}
		free(element);
		free(result);
		free(mul);
		fin.close();
		fout.close();

		return 0;
	}
	else {

	/* knn cuda */
	cudaError_t r;
	int *d_element, *d_result, *d_mul;
	r = cudaMalloc((void**)&d_element, sizeof(int)*MAXM*MAXN);
	// printf(" Malloc d_element : %s\n", cudaGetErrorString(r));
	r = cudaMalloc((void**)&d_result, sizeof(int)*MAXM*MAXK);
	// printf(" Malloc d_result : %s\n", cudaGetErrorString(r));
	r = cudaMalloc((void**)&d_mul, sizeof(int)*MAXM*MAXM);
	// printf(" Malloc d_mul : %s\n", cudaGetErrorString(r));

	/* time */
	cudaEvent_t begin, stop;
	cudaEventCreate(&begin);
	cudaEventCreate(&stop);
	cudaEventRecord(begin, 0);	//get begin time

	r = cudaMemcpy(d_element, element, sizeof(int)*MAXM*MAXN, cudaMemcpyHostToDevice);
	// printf(" Memory Copy d_element : %s\n", cudaGetErrorString(r));
	// r = cudaMemcpy(d_result, result, sizeof(int)*MAXM*MAXK, cudaMemcpyHostToDevice);
	// printf(" Memory Copy d_result : %s\n", cudaGetErrorString(r));
//
	// knn<<<BLOCKSPERGRID, THREADSPERBLOCK>>>(interger, &para);
	// getSum1<<<dim3(m/BLOCK, m/BLOCK, 1), dim3(BLOCK, BLOCK, 1)>>>(m, n, d_element, d_mul);
	getSum2<<<dim3(m/BLOCK, m/BLOCK, 1), dim3(BLOCK, BLOCK, 1)>>>(m, n, d_element, d_mul);
	cudaThreadSynchronize();
	getknn5<<<m/32, 32>>>(m, k, d_mul, d_result);

	r = cudaMemcpy(result, d_result, sizeof(int)*MAXM*MAXK, cudaMemcpyDeviceToHost);
	// printf(" Memcpy result : %s\n", cudaGetErrorString(r));

	cudaEventRecord(stop, 0);	//get stop time
	cudaEventSynchronize(stop);

	/* output */
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < k; j++) {
			printf("%d ", result[i*k+j]);
			fout << result[i*k+j] << ' ';
		}
		printf("\n");
		fout << endl;
	}

	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, begin, stop);
	// printf(" Excuetion Time on GPU: elapsedTime = %3.20f s\n", elapsedTime/1000);
	switch(m) {
		case 16384:
			printf("LARGE:<%3.20f s>\n", elapsedTime/1000);
			fout << "LARGE:<" << elapsedTime/1000  << " s>" << endl;
			break;
		case 4096:
			printf("MIDDLE:<%3.20f s>\n", elapsedTime/1000);
			fout << "MIDDLE:<" << elapsedTime/1000  << " s>" << endl;
			break;
		case 1024:
			printf("SMALL:<%3.20f s>\n", elapsedTime/1000);
			fout << "SMALL:<" << elapsedTime/1000  << " s>" << endl;
			break;
		default:
			break;
	}

	free(element);
	free(result);
	free(mul);
	cudaFree(d_element);
	cudaFree(d_result);
	cudaFree(d_mul);
	fin.close();
	fout.close();

	return 0;
}
}