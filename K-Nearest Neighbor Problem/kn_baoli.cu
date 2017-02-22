//kn.cu
//nvcc -arch=sm_30 kn.cu -o kn.out
/*
when m = 1024, time = 0.81 s
when m = 4096, time = 45.42 s
when m = 16384, time = 2875.6 s (getSum:2866.6s + getknn:2.78s)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>	//memset
#include <fstream>
// #include <iostream>
#include <time.h>
using namespace std;

#define MAXM 16384
#define MAXN 4096
#define MAXK 32
#define MAX_VALUE 5e5	//20^2*1000 = 4e5

int m, n, k, element[MAXM][MAXN], mul[MAXM][MAXM], result[MAXM][MAXK];

int getSum(int i, int j) {
	int tmpSum = 0;
	for (int g = 0; g < n; g++) {
		tmpSum += ((element[i][g] - element[j][g]) * (element[i][g] - element[j][g]));
	}
	return tmpSum;
}

void getknn(int tt) {
	int i = tt;
	int tmpMul[32];

	for (int j = 0; j < k; j++) {
		tmpMul[j] = MAX_VALUE;
	}

	for (int j = 0; j < m; j++) {
		int g;
		for(g = k-1; g >= 0; g--) {
			if(mul[i][j] >= tmpMul[g]) {
				break;
			}
		}
		if(g < k - 1) {
			for(int v = k - 2; v > g; v--) {
				result[i][v+1] = result[i][v];
				tmpMul[v+1] = tmpMul[v];
			}
			result[i][g+1] = j;
			tmpMul[g+1] = mul[i][j];
		}
	}
}

int main(int argc, char const *argv[])
{

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
	fin >> m >> n >> k;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			fin >> element[i][j];
	
printf("%d\n", __LINE__);
	clock_t start, stop, getSum1, getSum2, getknn1, getknn2;
    start = clock();
	/* main */
	getSum1 = clock();
	for (int i = 0; i < m; i++) {
		// printf("getSum i = %d\n", i);
		mul[i][i] = MAX_VALUE;
		for (int j = i+1; j < m; j++) {
			mul[i][j] = getSum(i, j);
			mul[j][i] = mul[i][j];
			// printf("mul[%d][%d] = %d\n", i, j, mul[i][j]);
		}
	}

/* example: testmul_128.out 
for (int i = 0; i < m; ++i)
{
	for(int j = 0; j < m; ++j) {
		fout << "mul[" << i << "][" << j << "] = " << mul[i][j] << endl;
	}
}*/

	getSum2 = clock();
// printf("%d\n", __LINE__);
	getknn1 = clock();
	for (int i = 0; i < m; i++) {
		// printf("getknn i = %d\n", i);
		getknn(i);
	}
	getknn2 = clock();
// printf("%d\n", __LINE__);
	/* output */
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < k; j++) {
			fout << result[i][j] << ' ';
		}
		fout << endl;
	}
printf("%d\n", __LINE__);
	fin.close();
	fout.close();

	stop = clock();
	printf("time getSum: %.5lf s\n", (double)(getSum2 - getSum1) / CLOCKS_PER_SEC);
	printf("time getknn: %.5lf s\n", (double)(getknn2 - getknn1) / CLOCKS_PER_SEC);
	printf("time elapsed: %.5lf s\n", (double)(stop - start) / CLOCKS_PER_SEC);

	return 0;
}