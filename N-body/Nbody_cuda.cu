//Nbody_cuda.cu
//g++ Nbody_pthread.c -lX11 -lrt -lpthread -o Nbody_pthread.out
//nvcc -arch=sm_20 Nbody_cuda.cu -lX11 -o Nbody_cuda.out

#include <X11/Xlib.h>	//-lX11
#include <stdio.h>
// #include <pthread.h>	//-lpthread
#include <stdio.h>		//freopen
#include <stdlib.h>
#include <time.h>	//-lrt
#include <math.h>
#include <string.h>	//memset
#include <cuda_runtime.h>
#include <iostream>
using namespace std;
#define MAXN 7300000	//max of N

// #define MAXN 1000000	//max of N
double G = 6.67e-11;  	//gravitational constant

int NUM_THREADS;//#NUM_THREADS: number of threads
double m;		//m: mass of each celestial body(float)
int T;			//T: number of steps
double t;		//t: time between each step (float)
//string FILE;	//FILE: your input file
int N;			//N: number of celestial bodies
double *x, *y;	//the initial position of each body
double *vx, *vy;	//the initial velocity of each body
double theta;	//θ: use in Barnes-Hut Algorithm
int enable;		//enable/disable: enable or disable Xwindow
double xmin, ymin;			//xmin, ymin: the lower left coordinate of Xwindow
double length;		//length: the length of the coordinate axis
double Length;		//Length: the Length of Window’s side (will be 10n times of length)
double Xmin, Xmax, Ymin, Ymax;

//body
double *Fx, *Fy;

//node
int *father;
int *child;
// int **child;
int *cnt;		//count the num of bodies
double *Sm;
double *Ax, *Ay;
double *X1, *X2;
double *Y1, *Y2;
int maxIndex = 1;
int Index = 0;
int childId;	//0, 1, 2, 3

//Window
Display *display;
Window window;
int screen;
GC gc;
XGCValues values;
int width, height;	//window (= Length)

void drawWindow() {
	printf(" drawWindow init\n");
	display = XOpenDisplay(NULL);
	if(display == NULL) {
		fprintf(stderr, "cannot open display\n");
		return;
	}

	int whiteColor = WhitePixel( display, DefaultScreen( display ) );
	int blackColor = BlackPixel( display, DefaultScreen( display ) );

	screen = DefaultScreen(display);

	/* set window size */
	width = Length;
	height = Length;

	/* set window position */
	int window_x = 0;
	int window_y = 0;

	/* border width in pixels */
	int border_width = 0;

	/* create window */
	window = XCreateSimpleWindow(display, RootWindow(display, screen), window_x, window_y, width, height, border_width, blackColor, whiteColor);

	/* create graph */
	long valuemask = 0;

	gc = XCreateGC(display, window, valuemask, &values);
	// XSetBackground (display, gc, whiteColor);
	XSetForeground (display, gc, blackColor);
	XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);

	/* map(show) the window */
	XMapWindow(display, window);
	XSync(display, 0);
}

bool inWindow(int i) {
	return (x[i] >= xmin && x[i] <= xmin + length && y[i] >= ymin && y[i] <= ymin + length);
}

void updateWindow() {
	XClearWindow(display, window);
	int i;
	for(i = 0; i < N; i++) {
		if (inWindow(i)) {
			XFillRectangle(display, window, gc, (x[i]-xmin)/length*Length-1, (y[i]-ymin)/length*Length-1, 2, 2);
		}
	}
	// XFlush(display);
	XSync(display,true);
}

// int g = 4;
void putBody(int i) {
	while (1){
		if (X2[Index] - X1[Index] < 1e-10) {
			cnt[Index]++;
			Ax[Index] = (Ax[Index]*Sm[Index]+x[i] * m) / (Sm[Index]+m);
			Ay[Index] = (Ay[Index]*Sm[Index]+y[i] * m) / (Sm[Index]+m);
			Sm[Index] += m;
			break;
		}

		//1. If node x does not contain a body, put the new body b here.
		if (cnt[Index] == 0) {
			cnt[Index]++;
			Ax[Index] = x[i];
			Ay[Index] = y[i];
			Sm[Index] = m;
			break;
		}

		//2. If node x is an internal node, update the center-of-mass and total mass of x.
		//Recursively insert the body b in the appropriate quadrant.
		if (cnt[Index] > 1) {
			cnt[Index]++;
			Ax[Index] = (Ax[Index] * Sm[Index] + x[i] * m) / (Sm[Index] + m);
			Ay[Index] = (Ay[Index] * Sm[Index] + y[i] * m) / (Sm[Index] + m);
			Sm[Index] += m;
			if (2*x[i] <= X1[Index]+X2[Index])	childId = 0;
			else	childId = 2;
			if (2*y[i] > Y1[Index]+Y2[Index])	childId += 1;
			// Index = child[Index][childId];
			Index = child[4*Index+childId];
			continue;
		}

		//3. If node x is an external node, say containing a body named c, then there are two bodies b and c in the same region.
		//Subdivide the region further by creating four children.
		//Then, recursively insert both b and c into the appropriate quadrant(s).
		//Since b and c may still end up in the same quadrant, there may be several subdivisions during a single insertion.
		//Finally, update the center-of-mass and total mass of x.
		if (cnt[Index] == 1) {
			cnt[Index]++;

			//create four children
			child[4*Index] = maxIndex;
			// child[Index][0] = maxIndex;
			father[maxIndex] = Index;
			X1[maxIndex] = X1[Index];
			X2[maxIndex] = (X1[Index] + X2[Index]) / 2;
			Y1[maxIndex] = Y1[Index];
			Y2[maxIndex] = (Y1[Index] + Y2[Index]) / 2;
			maxIndex++;

			child[4*Index+1] = maxIndex;
			// child[Index][1] = maxIndex;
			father[maxIndex] = Index;
			X1[maxIndex] = X1[Index];
			X2[maxIndex] = (X1[Index] + X2[Index]) / 2;
			Y1[maxIndex] = (Y1[Index] + Y2[Index]) / 2;
			Y2[maxIndex] = Y2[Index];
			maxIndex++;

			child[4*Index+2] = maxIndex;
			// child[Index][2] = maxIndex;
			father[maxIndex] = Index;
			X1[maxIndex] = (X1[Index] + X2[Index]) / 2;
			X2[maxIndex] = X2[Index];
			Y1[maxIndex] = Y1[Index];
			Y2[maxIndex] = (Y1[Index] + Y2[Index]) / 2;
			maxIndex++;

			child[4*Index+3] = maxIndex;
			// child[Index][3] = maxIndex;
			father[maxIndex] = Index;
			X1[maxIndex] = (X1[Index] + X2[Index]) / 2;
			X2[maxIndex] = X2[Index];
			Y1[maxIndex] = (Y1[Index] + Y2[Index]) / 2;
			Y2[maxIndex] = Y2[Index];
			maxIndex++;
			for(int j = 0; j < 4; j++) {
				// if(child[Index][j] < Index)
				if(child[4*Index+j] < Index)
				{
					printf("Erorr\n");
					printf(" e Index = %d\n", Index);
					printf(" e maxIndex = %d\n", maxIndex);
					sleep(100);
				}
			}

			//for origin body c
			if (2*Ax[Index] <= X1[Index]+X2[Index])	childId = 0;
			else	childId = 2;
			if (2*Ay[Index] > Y1[Index]+Y2[Index])	childId += 1;
			// cnt[child[Index][childId]]++;
			// Ax[child[Index][childId]] = Ax[Index];
			// Ay[child[Index][childId]] = Ay[Index];
			// Sm[child[Index][childId]] = Sm[Index];
			cnt[child[Index*4+childId]]++;
			Ax[child[Index*4+childId]] = Ax[Index];
			Ay[child[Index*4+childId]] = Ay[Index];
			Sm[child[Index*4+childId]] = Sm[Index];

			//for origin Index
			Ax[Index] = (Ax[Index] * Sm[Index] + x[i] * m) / (Sm[Index] + m);
			Ay[Index] = (Ay[Index] * Sm[Index] + y[i] * m) / (Sm[Index] + m);
			Sm[Index] += m;

			//for new body b
			if (2*x[i] <= X1[Index]+X2[Index])	childId = 0;
			else	childId = 2;
			if (2*y[i] > Y1[Index]+Y2[Index])	childId += 1;
			// Index = child[Index][childId];
			Index = child[Index*4+childId];
		}
	}
}

__global__ void calc(double m, double G, double theta, double* x, double* y, double* Fx, double* Fy, int* child, int* cnt, double* Sm, double* Ax, double* Ay, double* X1, double* X2) {
	printf("###### gpu***\n");
	
	int stack[100000000];//...
	
	// printf("\n\ncnt[0] = %d\n", cnt[0]);
	// printf("x[0] = %f\n", x[0]);
	// printf("y[0] = %f\n", y[0]);
	// printf("Fx[0] = %f\n", Fx[0]);
	// printf("Fy[0] = %f\n", Fy[0]);
	// printf("Sm[0] = %f\n", Sm[0]);
	// printf("Ax[0] = %f\n", Ax[0]);
	// printf("Ay[0] = %f\n", Ay[0]);
	// printf("X1[0] = %f\n", X1[0]);

	// int *stack;
	// cudaError_t r;
	// r = cudaMalloc((void**)&stack, sizeof(int)*5*MAXN);
	// printf(" Malloc stack : %s\n", cudaGetErrorString(r));

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
		//calcForce
		int pi = 1, pIndex;
		stack[1] = 0;
		// pi = 9;//
		while(pi>0) {
			pIndex = stack[pi];
			
			if (cnt[pIndex] == 0) {
				// printf("  0\n");
				pi--;
				continue;
			}
			
			if (cnt[pIndex] == 1) {
				
				// printf("  1\n");
				double r2 = (Ax[pIndex]-x[i])*(Ax[pIndex]-x[i])+(Ay[pIndex]-y[i])*(Ay[pIndex]-y[i]);
				double r = sqrt(r2);
				if(r < 1e-10) {pi--; continue;}
				double r3 = r2*r;
				Fx[i] += G*m*Sm[pIndex]*(Ax[pIndex]-x[i])/r3;
				Fy[i] += G*m*Sm[pIndex]*(Ay[pIndex]-y[i])/r3;
				pi--;
				continue;
			}

			if (cnt[pIndex] > 1) {
				// printf("  2\n");
				//calc s/d
				double d2 = (Ax[pIndex]-x[i])*(Ax[pIndex]-x[i])+(Ay[pIndex]-y[i])*(Ay[pIndex]-y[i]);
				double d = sqrt(d2);
				if(d < 1e-10)	{pi--; continue;}
				double d3 = d2*d;
				
				// if ((X2[pIndex]-X1[pIndex] < d * theta) || (X2[pIndex] - X1[pIndex] < 1e-10)) {
				// if (1) {
					Fx[i] += G*m*Sm[pIndex]*(Ax[pIndex]-x[i])/d3;
					Fy[i] += G*m*Sm[pIndex]*(Ay[pIndex]-y[i])/d3;
					pi--;
					continue;
				// }
				// else {
				// 	Fx[i] += G*m*Sm[pIndex]*(Ax[pIndex]-x[i])/d3;
				// 	Fy[i] += G*m*Sm[pIndex]*(Ay[pIndex]-y[i])/d3;
				// 	pi--;
				// 	continue;
					// printf("pi = %d\n", pi);
					// stack[pi+3] = child[4*pIndex+3];
					// stack[pi+2] = child[4*pIndex+2];
					// stack[pi+1] = child[4*pIndex+1];
					// stack[pi] = child[4*pIndex];
					// pi += 3;
					// continue;
				// }

			}

		}
		// cudaFree(stack);
	// free stack;
		// printf("####### gpu i = %d, Fx = %f, Fy = %f\n", i, Fx[i], Fy[i]);
	
}

int main(int argc, char const *argv[])
{
	cudaSetDevice(1);
	cudaError_t r;
	if (!(argc == 12 || argc == 8)) {
		printf(" Input: ./a.out #threads m T t FILE theta enable/disable xmin ymin length Length\n");
		printf(" Example:\n");
		printf(" test1.txt (N=945): ./a.out 12 1 1000 0.5 test1.txt 0 enable -0.3 -0.3 0.6 600\n");
		printf(" test2.txt (N=81921): ./a.out 12 1 1000 0.01 test2.txt 1 enable 0 -20 1 500\n");
		printf(" test3.txt (N=1923840): ./a.out 12 1 1000 0.01 test3.txt 1 enable 0 0 1 500\n");
		printf(" test4.txt (N=7283942): ./a.out 12 1 1000 0.001 test4.txt 1 enable 0 0 1 500\n");
		return 0;
	}
	//malloc
	int i;
	x = (double*)malloc(sizeof(double)*MAXN);
	y = (double*)malloc(sizeof(double)*MAXN);
	vx = (double*)malloc(sizeof(double)*MAXN);
	vy = (double*)malloc(sizeof(double)*MAXN);
	Fx = (double*)malloc(sizeof(double)*MAXN);
	Fy = (double*)malloc(sizeof(double)*MAXN);

	father = (int*)malloc(sizeof(int)*5*MAXN);
	// child = (int**)malloc(sizeof(int *)*5*MAXN);
	// for (i = 0; i < 5*MAXN; i++)
	// 	child[i] = (int*)malloc(sizeof(int)*4);
	child = (int*)malloc(sizeof(int)*(5*MAXN+10)*4);
	cnt = (int*)malloc(sizeof(int)*5*MAXN);
	Sm = (double*)malloc(sizeof(double)*5*MAXN);
	Ax = (double*)malloc(sizeof(double)*5*MAXN);
	Ay = (double*)malloc(sizeof(double)*5*MAXN);
	X1 = (double*)malloc(sizeof(double)*5*MAXN);
	X2 = (double*)malloc(sizeof(double)*5*MAXN);
	Y1 = (double*)malloc(sizeof(double)*5*MAXN);
	Y2 = (double*)malloc(sizeof(double)*5*MAXN);
printf(" 1. malloc over\n");

		//time
		struct timespec time11 = {0, 0};
		struct timespec time12 = {0, 0};
		struct timespec time21 = {0, 0};
		struct timespec time22 = {0, 0};
		// struct timespec time31 = {0, 0};
		// struct timespec time32 = {0, 0};
		struct timespec time41 = {0, 0};
		struct timespec time42 = {0, 0};
		clock_gettime(CLOCK_MONOTONIC, &time11);

	//ReadIn
	NUM_THREADS = atoi(argv[1]);
	m = atof(argv[2]);
	T = atoi(argv[3]);
	t = atof(argv[4]);
	freopen(argv[5], "r", stdin);	//FILE
	//freopen("output.txt", "w", stdout);
	cin >> N;
	int once = 1;
	for (i = 0; i < N; ++i) {
		cin >> x[i] >> y[i] >> vx[i] >> vy[i];
		if (once) {
			Xmin = Xmax = x[i];
			Ymin = Ymax = y[i];
			once = 0;
		}
		if (x[i] < Xmin)	Xmin = x[i];
		else if (x[i] > Xmax)	Xmax = x[i];
		if (y[i] < Ymin)	Ymin = y[i];
		else if (y[i] > Ymax)	Ymax = y[i];
	}
	double dx = Xmax - Xmin;
	double dy = Ymax - Ymin;
	Xmin -= 0.1 * dx;
	Xmax += 0.1 * dx;
	Ymin -= 0.1 * dy;
	Ymax += 0.1 * dy;
	X1[0] = Xmin;
	X2[0] = Xmax;
	Y1[0] = Ymin;
	Y2[0] = Ymax;

	if (strcmp(argv[7], "enable") == 0)	enable = 1;
	else	enable = 0;
	if (enable)	{
		theta = atof(argv[6]);
		xmin = atof(argv[8]);
		ymin = atof(argv[9]);
		length = atof(argv[10]);
		Length = atof(argv[11]);
		//init window
		drawWindow();
	}
	clock_gettime(CLOCK_MONOTONIC, &time12);
printf(" 2. ReadIn over\n");

		//cudaMalloc
		double *d_x;
		r = cudaMalloc((void**)&d_x, sizeof(double)*MAXN);
		printf(" Malloc d_x : %s\n", cudaGetErrorString(r));

		double *d_y;
		r = cudaMalloc((void**)&d_y, sizeof(double)*MAXN);
		printf(" Malloc d_y : %s\n", cudaGetErrorString(r));

		double *d_Fx;
		r = cudaMalloc((void**)&d_Fx, sizeof(double)*MAXN);
		printf(" Malloc d_Fx : %s\n", cudaGetErrorString(r));

		double *d_Fy;
		r = cudaMalloc((void**)&d_Fy, sizeof(double)*MAXN);
		printf(" Malloc d_Fy : %s\n", cudaGetErrorString(r));

		int *d_child;
		r = cudaMalloc((void**)&d_child, sizeof(int)*5*MAXN*4);
		printf(" Malloc d_child : %s\n", cudaGetErrorString(r));

		int *d_cnt;
		r = cudaMalloc((void**)&d_cnt, sizeof(int)*5*MAXN);
		printf(" Malloc d_cnt : %s\n", cudaGetErrorString(r));

		double *d_Sm;
		r = cudaMalloc((void**)&d_Sm, sizeof(double)*5*MAXN);
		printf(" Malloc d_Sm : %s\n", cudaGetErrorString(r));

		double *d_Ax;
		r = cudaMalloc((void**)&d_Ax, sizeof(double)*5*MAXN);
		printf(" Malloc d_Ax : %s\n", cudaGetErrorString(r));

		double *d_Ay;
		r = cudaMalloc((void**)&d_Ay, sizeof(double)*5*MAXN);
		printf(" Malloc d_Ay : %s\n", cudaGetErrorString(r));

		double *d_X1;
		r = cudaMalloc((void**)&d_X1, sizeof(double)*5*MAXN);
		printf(" Malloc d_X1 : %s\n", cudaGetErrorString(r));

		double *d_X2;
		r = cudaMalloc((void**)&d_X2, sizeof(double)*5*MAXN);
		printf(" Malloc d_X2 : %s\n", cudaGetErrorString(r));
printf(" 3. while T start...\n");


	// int *d_stack;
	// r = cudaMalloc((void**)&d_stack, sizeof(int)*5*MAXN);
	// printf(" Malloc d_stack : %s\n", cudaGetErrorString(r));


	//main
	while(T--) {
		printf(" *** T = %d ***\n", T);
		//init
		for (i = 0; i < N; i++) {
			Fx[i] = Fy[i] = 0;
		}
		for (i = 0; i < 5*MAXN; i++) {
			cnt[i] = Sm[i] = Ax[i] = Ay[i] = 0;
		}

		maxIndex = 1;

		clock_gettime(CLOCK_MONOTONIC, &time21);

		//buildTree
		//printf(" buildTree\n");
		for (i = 0; i < N; i++)	{
// printf(" ... T = %d, i = %d\n", T, i);
			Index = 0;
			putBody(i);
		}
		clock_gettime(CLOCK_MONOTONIC, &time22);
// printf("buildTree over, T = %d\n", T);

		//gpu
		int THREADSPERBLOCK = 512;
		int BLOCKSPERGRID = (int)ceil((double)N / 512);
		// int BLOCKSPERGRID = 1024;

		//calcForce

	r = cudaMemcpy(d_x, x, sizeof(double)*MAXN, cudaMemcpyHostToDevice);
	printf(" Memory Copy d_x : %s\n", cudaGetErrorString(r));
	r = cudaMemcpy(d_y, y, sizeof(double)*MAXN, cudaMemcpyHostToDevice);
	printf(" Memory Copy d_y : %s\n", cudaGetErrorString(r));
	r = cudaMemcpy(d_Fx, Fx, sizeof(double)*MAXN, cudaMemcpyHostToDevice);
	printf(" Memory Copy d_Fx : %s\n", cudaGetErrorString(r));
	r = cudaMemcpy(d_Fy, Fy, sizeof(double)*MAXN, cudaMemcpyHostToDevice);
	printf(" Memory Copy d_Fy : %s\n", cudaGetErrorString(r));
	r = cudaMemcpy(d_child, child, sizeof(int)*5*MAXN*4, cudaMemcpyHostToDevice);
	printf(" Memory Copy d_child : %s\n", cudaGetErrorString(r));
	r = cudaMemcpy(d_cnt, cnt, sizeof(int)*5*MAXN, cudaMemcpyHostToDevice);
	printf(" Memory Copy d_cnt : %s\n", cudaGetErrorString(r));
	r = cudaMemcpy(d_Sm, Sm, sizeof(double)*5*MAXN, cudaMemcpyHostToDevice);
	printf(" Memory Copy d_Sm : %s\n", cudaGetErrorString(r));
	r = cudaMemcpy(d_Ax, Ax, sizeof(double)*5*MAXN, cudaMemcpyHostToDevice);
	printf(" Memory Copy d_Ax : %s\n", cudaGetErrorString(r));
	r = cudaMemcpy(d_Ay, Ay, sizeof(double)*5*MAXN, cudaMemcpyHostToDevice);
	printf(" Memory Copy d_Ay : %s\n", cudaGetErrorString(r));
	r = cudaMemcpy(d_X1, X1, sizeof(double)*5*MAXN, cudaMemcpyHostToDevice);
	printf(" Memory Copy d_X1 : %s\n", cudaGetErrorString(r));
	r = cudaMemcpy(d_X2, X2, sizeof(double)*5*MAXN, cudaMemcpyHostToDevice);
	printf(" Memory Copy d_X2 : %s\n", cudaGetErrorString(r));
// printf("child\n");
		// printf("(o)d_x address = %ld\n", (long int)d_x);
		// printf("(o)x address = %ld\n", (long int)x);
// printf("(o)Fx[0] = %f\n", Fx[0]);//...
		//time
		cudaEvent_t begin, stop;
		cudaEventCreate(&begin);
		cudaEventCreate(&stop);	
		//Get begin time
		cudaEventRecord(begin, 0);
		// clock_gettime(CLOCK_MONOTONIC, &time2);
		calc<<<BLOCKSPERGRID, THREADSPERBLOCK>>>(m, G, theta, d_x, d_y, d_Fx, d_Fy, d_child, d_cnt, d_Sm, d_Ax, d_Ay, d_X1, d_X2);
		cudaThreadSynchronize();
		// clock_gettime(CLOCK_MONOTONIC, &time3);
	//Get stop time
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, begin, stop);
	printf(" CLOCK_calculate_TIME(Execution Time on GPU): %3.20f s\n", elapsedTime/1000);

	r = cudaMemcpy(Fx, d_Fx, sizeof(double)*MAXN, cudaMemcpyDeviceToHost);
	printf(" Memcpy Fx : %s\n", cudaGetErrorString(r));
	r = cudaMemcpy(Fy, d_Fy, sizeof(double)*MAXN, cudaMemcpyDeviceToHost);
	printf(" Memcpy Fy : %s\n", cudaGetErrorString(r));
	// r = cudaMemcpy(cnt, d_cnt, sizeof(int)*5*MAXN, cudaMemcpyDeviceToHost);
	// printf(" Memcpy cnt : %s\n", cudaGetErrorString(r));
// printf("(n)Fx[0] = %f\n", Fx[0]);//...

// printf("(n)Fx address = %ld\n", (long int)Fx);
// printf("(n)d_Fx address = %ld\n", (long int)d_Fx);


		//get new vx, vy, x, y
		for (i = 0; i < N; i++) {
			// printf("Fx[%d] = %f\n", i, Fx[i]);
			// printf("Fy[%d] = %f\n", i, Fy[i]);
			// sleep(1);
			vx[i] += Fx[i]*t/m;
			vy[i] += Fy[i]*t/m;
			x[i] += vx[i]*t;
			y[i] += vy[i]*t;
			if (x[i] < Xmin)	x[i] = 2*Xmin - x[i];
			if (x[i] > Xmax)	x[i] = 2*Xmax - x[i];
			if (y[i] < Ymin)	y[i] = 2*Ymin - y[i];
			if (y[i] > Ymax)	y[i] = 2*Ymax - y[i];
		}
		clock_gettime(CLOCK_MONOTONIC, &time41);
// printf(" joined!\n");
		if (enable) {
			updateWindow();
		}
		clock_gettime(CLOCK_MONOTONIC, &time42);
		printf(" CLOCK_IO_TIME: %lf s\n", (time12.tv_sec - time11.tv_sec) + (time12.tv_nsec - time11.tv_nsec) / 1e9);
		printf(" CLOCK_buildTree_TIME: %lf s\n", (time22.tv_sec - time21.tv_sec) + (time22.tv_nsec - time21.tv_nsec) / 1e9);
		// printf(" CLOCK_calculate_TIME: %lf s\n", (time32.tv_sec - time31.tv_sec) + (time32.tv_nsec - time31.tv_nsec) / 1e9);
		printf(" CLOCK_drawWindow_TIME: %lf s\n", (time42.tv_sec - time41.tv_sec) + (time42.tv_nsec - time41.tv_nsec) / 1e9);

	}

    //free
    free(x);
    free(y);
    free(vx);
    free(vy);
    free(Fx);
    free(Fy);
    free(father);
  //   for (i = 0; i < 5*MAXN; i++)
		// free(child[i]);
	free(child);
	free(cnt);
	free(Sm);
	free(Ax);
	free(Ay);
	free(X1);
	free(X2);
	free(Y1);
	free(Y2);

	cudaFree(d_x);
	cudaFree(d_y);
	cudaFree(d_Fx);
	cudaFree(d_Fy);
	cudaFree(d_child);
	cudaFree(d_cnt);
	cudaFree(d_Sm);
	cudaFree(d_Ax);
	cudaFree(d_Ay);
	cudaFree(d_X1);
	cudaFree(d_X2);

	return 0;
}