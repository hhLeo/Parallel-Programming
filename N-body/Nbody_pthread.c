//Nbody_pthread.c
//gcc Nbody_pthread.c -lX11 -lrt -lpthread -o Nbody_pthread.out
//time =  ms

#include <X11/Xlib.h>	//-lX11
#include <stdio.h>
#include <pthread.h>	//-lpthread
#include <stdio.h>		//freopen
#include <stdlib.h>
#include <time.h>	//-lrt
#include <math.h>
#include <string.h>	//memset
#include <iostream>
using namespace std;

#define MAXN 10000000	//max of N
#define G 6.67e-11  			//gravitational constant

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
int **child;
int *cnt;		//count the num of bodies
double *Sm;
double *Ax, *Ay;
double *X1, *X2;
double *Y1, *Y2;
int maxIndex = 1;
int Index = 0;
int childId;	//0, 1, 2, 3

//pthread
pthread_mutex_t count_mutex;
int nowN;

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
	// printf("update...\n");
	// XFillRectangle(display, window, gc, 20, 20, 10, 10);
	// XFillArc(display, window, gc, 250, 250, 10, 10, 0, 360);
	int i, j;
	for(i = 0; i < N; i++) {
		if (inWindow(i)) {
			// printf("x[i] = %f, y[i] = %f\n", x[i], y[i]);
			// XFillRectangle(display, window, gc, 20, 20, 5, 5);
			//XFillArc(display, d, gc, x, y, width, height, angle1, angle2);
			XFillRectangle(display, window, gc, (x[i]-xmin)/length*Length-1, (y[i]-ymin)/length*Length-1, 2, 2);
		}
	}
	// XFlush(display);
	XSync(display,true);
	// sleep(5);
}

// int g = 4;
void putBody(int i) {
	// printf("i = %d, x = %f, y = %f\n", i, x[i], y[i]);
	while (1){
		// Ax[Index] = (Ax[Index]*Sm[Index]+x[i] * m) / (Sm[Index]+m);
		// Ay[Index] = (Ay[Index]*Sm[Index]+y[i] * m) / (Sm[Index]+m);
		// Sm[Index] += m;
		if (X2[Index] - X1[Index] < 1e-10) {
			cnt[Index]++;
			Ax[Index] = (Ax[Index]*Sm[Index]+x[i] * m) / (Sm[Index]+m);
			Ay[Index] = (Ay[Index]*Sm[Index]+y[i] * m) / (Sm[Index]+m);
			Sm[Index] += m;
			break;
		}

		//1. If node x does not contain a body, put the new body b here.
		if (cnt[Index] == 0) {
// printf("cnt=0; cnt[%d]++\n", Index);
// printf("%d\n", __LINE__);
			cnt[Index]++;
			// Ax[Index] = (Ax[Index] * Sm[Index] + x[i] * m) / (Sm[Index] + m);
			// Ay[Index] = (Ay[Index] * Sm[Index] + y[i] * m) / (Sm[Index] + m);
			// Sm[Index] += m;
			Ax[Index] = x[i];
			Ay[Index] = y[i];
			Sm[Index] = m;
			break;
		}

		//2. If node x is an internal node, update the center-of-mass and total mass of x.
		//Recursively insert the body b in the appropriate quadrant.
		if (cnt[Index] > 1) {
// printf("cnt[%d] = %d\n", Index, cnt[Index]);
			// if (g--)
			// 	printf("cnt>1; cnt[%d]++\n", Index);
// printf("%d\n", __LINE__);
			cnt[Index]++;
			Ax[Index] = (Ax[Index] * Sm[Index] + x[i] * m) / (Sm[Index] + m);
			Ay[Index] = (Ay[Index] * Sm[Index] + y[i] * m) / (Sm[Index] + m);
			Sm[Index] += m;
			if (2*x[i] <= X1[Index]+X2[Index])	childId = 0;
			else	childId = 2;
			if (2*y[i] <= Y1[Index]+Y2[Index])	childId += 0;
			else	childId += 1;
// printf("(orig)Index = %d\n", Index);
			Index = child[Index][childId];
// printf("(new)Index = %d\n", Index);
			continue;
		}

		//3. If node x is an external node, say containing a body named c, then there are two bodies b and c in the same region.
		//Subdivide the region further by creating four children.
		//Then, recursively insert both b and c into the appropriate quadrant(s).
		//Since b and c may still end up in the same quadrant, there may be several subdivisions during a single insertion.
		//Finally, update the center-of-mass and total mass of x.
		if (cnt[Index] == 1) {
			// if (x[i] == Ax[Index] && y[i] == Ay[Index])	x[i] += 1e3;
// printf("%d\n", __LINE__);
// printf("x = %f, Ax = %f, y = %f, Ay = %f\n", x[i], Ax[Index], y[i], Ay[Index]);
// printf(" s Index = %d\n", Index);
// printf(" s maxIndex = %d\n", maxIndex);

// printf("cnt=1; cnt[%d]++\n", Index);
			cnt[Index]++;

			//create four children
			child[Index][0] = maxIndex;
			father[maxIndex] = Index;
			X1[maxIndex] = X1[Index];
			X2[maxIndex] = (X1[Index] + X2[Index]) / 2;
			Y1[maxIndex] = Y1[Index];
			Y2[maxIndex] = (Y1[Index] + Y2[Index]) / 2;
			maxIndex++;

			child[Index][1] = maxIndex;
			father[maxIndex] = Index;
			X1[maxIndex] = X1[Index];
			X2[maxIndex] = (X1[Index] + X2[Index]) / 2;
			Y1[maxIndex] = (Y1[Index] + Y2[Index]) / 2;
			Y2[maxIndex] = Y2[Index];
			maxIndex++;

			child[Index][2] = maxIndex;
			father[maxIndex] = Index;
			X1[maxIndex] = (X1[Index] + X2[Index]) / 2;
			X2[maxIndex] = X2[Index];
			Y1[maxIndex] = Y1[Index];
			Y2[maxIndex] = (Y1[Index] + Y2[Index]) / 2;
			maxIndex++;

			child[Index][3] = maxIndex;
			father[maxIndex] = Index;
			X1[maxIndex] = (X1[Index] + X2[Index]) / 2;
			X2[maxIndex] = X2[Index];
			Y1[maxIndex] = (Y1[Index] + Y2[Index]) / 2;
			Y2[maxIndex] = Y2[Index];
			maxIndex++;
			for(int j = 0; j < 4; j++) {
				if(child[Index][j] < Index)
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
			if (2*Ay[Index] <= Y1[Index]+Y2[Index])	childId += 0;
			else	childId += 1;
// printf("cnt=1; (orig)cnt[%d]++\n", child[Index][childId]);
			cnt[child[Index][childId]]++;
			Ax[child[Index][childId]] = Ax[Index];
			Ay[child[Index][childId]] = Ay[Index];
			Sm[child[Index][childId]] = Sm[Index];

			//for origin Index
			Ax[Index] = (Ax[Index] * Sm[Index] + x[i] * m) / (Sm[Index] + m);
			Ay[Index] = (Ay[Index] * Sm[Index] + y[i] * m) / (Sm[Index] + m);
			Sm[Index] += m;

			//for new body b
			if (2*x[i] <= X1[Index]+X2[Index])	childId = 0;
			else	childId = 2;
			// if (2*y[i] <= Y1[Index]+Y2[Index])	childId += 0;
			// else	childId += 1;
			if (2*y[i] > Y1[Index]+Y2[Index])	childId += 1;
			Index = child[Index][childId];
// printf("nextIndex = %d\n", Index);
		}
	}
}

/*
//orig
void calcForce(int i, int pIndex) {
	if (cnt[pIndex] == 0)	return;
	if (cnt[pIndex] == 1) {
		double r2 = (Ax[pIndex]-x[i])*(Ax[pIndex]-x[i])+(Ay[pIndex]-y[i])*(Ay[pIndex]-y[i]);
		double r = sqrt(r2);
		if(r < 1e-10)	return;
		double r3 = r2*r;
		Fx[i] += G*m*m*(Ax[pIndex]-x[i])/r3;
		Fy[i] += G*m*m*(Ay[pIndex]-y[i])/r3;
		return;
	}
	if (cnt[pIndex] > 1) {
		//calc s/d
		double d2 = (Ax[pIndex]-x[i])*(Ax[pIndex]-x[i])+(Ay[pIndex]-y[i])*(Ay[pIndex]-y[i]);
		double d = sqrt(d2);
		if(d < 1e-10)	return;
		double d3 = d2*d;
		if ((X2[pIndex]-X1[pIndex])/sqrt(d2) < theta) {
			Fx[i] += G*m*m*(Ax[pIndex]-x[i])/d3;
			Fy[i] += G*m*m*(Ay[pIndex]-y[i])/d3;
			return;
		}
		else {
			int j;
			for(j = 0; j < 4; j++)
				calcForce(i, child[pIndex][j]);
			return;
		}
	}
}
*/

void *calc(void *t1) {
	int i;
	int *stack;
	stack = (int*)malloc(sizeof(int)*5*MAXN);
	while (nowN < N) {
		pthread_mutex_lock(&count_mutex);
		i = nowN++;
		pthread_mutex_unlock(&count_mutex);

		//calcForce(i, 0);	//orig
		int pi = 1;
		stack[1] = 0;
		while(pi>0) {
			if (cnt[stack[pi]] == 0) {
				pi--;
				continue;
			}
			if (cnt[stack[pi]] == 1) {
				double r2 = (Ax[stack[pi]]-x[i])*(Ax[stack[pi]]-x[i])+(Ay[stack[pi]]-y[i])*(Ay[stack[pi]]-y[i]);
				double r = sqrt(r2);
				if(r < 1e-10) {pi--; continue;}
				double r3 = r2*r;
				Fx[i] += G*m*Sm[stack[pi]]*(Ax[stack[pi]]-x[i])/r3;
				Fy[i] += G*m*Sm[stack[pi]]*(Ay[stack[pi]]-y[i])/r3;
				pi--;
				continue;
			}
			if (cnt[stack[pi]] > 1) {
				//calc s/d
				double d2 = (Ax[stack[pi]]-x[i])*(Ax[stack[pi]]-x[i])+(Ay[stack[pi]]-y[i])*(Ay[stack[pi]]-y[i]);
				double d = sqrt(d2);
				if(d < 1e-10)	{pi--; continue;}
				double d3 = d2*d;
				if ((X2[stack[pi]]-X1[stack[pi]]) < d * theta || X2[stack[pi]] - X1[stack[pi]] < 1e-10) {
					Fx[i] += G*m*Sm[stack[pi]]*(Ax[stack[pi]]-x[i])/d3;
					Fy[i] += G*m*Sm[stack[pi]]*(Ay[stack[pi]]-y[i])/d3;
					pi--;
					continue;
				}
				else {
					stack[pi+3] = child[stack[pi]][3];
					stack[pi+2] = child[stack[pi]][2];
					stack[pi+1] = child[stack[pi]][1];
					stack[pi] = child[stack[pi]][0];
					pi += 3;
					continue;
				}
			}
		}
	}
	// free stack;
	pthread_exit(NULL);
}

int main(int argc, char const *argv[])
{
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
	int i, j;
	x = (double*)malloc(sizeof(double)*MAXN);
	y = (double*)malloc(sizeof(double)*MAXN);
	vx = (double*)malloc(sizeof(double)*MAXN);
	vy = (double*)malloc(sizeof(double)*MAXN);
	Fx = (double*)malloc(sizeof(double)*MAXN);
	Fy = (double*)malloc(sizeof(double)*MAXN);

	father = (int*)malloc(sizeof(int)*5*MAXN);
	child = (int**)malloc(sizeof(int *)*5*MAXN);
	for (i = 0; i < 5*MAXN; i++)
		child[i] = (int*)malloc(sizeof(int)*4);
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
		struct timespec time31 = {0, 0};
		struct timespec time32 = {0, 0};
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
// printf("%s\n", argv[7]);

	if (strcmp(argv[7], "enable") == 0)	enable = 1;
	else	enable = 0;
// printf("enable = %d\n", enable);
	// sleep(3);
	if (enable)	{
		//readIn
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

	//threads
	pthread_t threads[NUM_THREADS];
// printf(" 3. while T start...\n");
	//main
	// T = 3;
	while(T--) {
		printf(" *** T = %d ***\n", T);
		if (T % 100 == 0)
			sleep(10);
		//init
		for (i = 0; i < N; i++) {
			Fx[i] = Fy[i] = 0;
		}
// printf("338\n");
		// for (i = 1; i < 5*MAXN; i++) {
		// 	X1[i] = X2[i] = Y1[i] = Y2[i] = 0;
		// }
// printf("342\n");
		// printf("memset\n");
		// memset(Sm, 0, sizeof(double)*5*MAXN);
		// memset(Ax, 0, sizeof(double)*5*MAXN);
		// memset(Ay, 0, sizeof(double)*5*MAXN);
		for (i = 0; i < 5*MAXN; i++) {
			cnt[i] = Sm[i] = Ax[i] = Ay[i] = 0;
		}
// printf("350\n");

		nowN = 0;	//for dynamic pthreads
		maxIndex = 1;
		// printf(" buildTree\n");

		clock_gettime(CLOCK_MONOTONIC, &time21);
		//buildTree
		for (i = 0; i < N; i++)	{
// printf(" ... T = %d, i = %d\n", T, i);
			Index = 0;
			putBody(i);
// printf(" maxIndex = %d\n", maxIndex);
			// for (j = 0; j < maxIndex; j++) {
			// 	printf("  cnt[%d] = %d\n", j, cnt[j]);
			// }
		}
		clock_gettime(CLOCK_MONOTONIC, &time22);
		
// printf("buildTree over, T = %d\n", T);
		// for (i = 0; i < N; ++i)
		// {
		// 	if(int(Sm[i]) != cnt[i])
		// 		printf("\n\nSm != cnt !!!! Sm[%d] = %f, cnt = %d\n\n", i, Sm[i], cnt[i]);
		// }
		// for(i = 0; i < maxIndex; i++) {
//		// 	printf("cnt[%d] = %d, ..\n", i, cnt[i]);
		// }
		// sleep(100);
		// for (i = 0; i < 10; ++i)
		// {
		// 	printf("x[%d] = %f, y[%d] = %f\n", i, x[i], i, y[i]);
		// }
		// sleep(10);
		clock_gettime(CLOCK_MONOTONIC, &time31);
		for (i = 0; i < NUM_THREADS; i++) {
			pthread_create(&threads[i], NULL, calc, (void *)((long int)i));
		}

		for (i = 0; i < NUM_THREADS; i++) {
	        pthread_join(threads[i], NULL);
		}

		// printf("t is %f\n", t);
		for (i = 0; i < N; i++) {
			vx[i] += Fx[i]*t/m;
			vy[i] += Fy[i]*t/m;
			x[i] += vx[i]*t;
			y[i] += vy[i]*t;
			if (x[i] < Xmin)	x[i] = 2*Xmin - x[i];
			if (x[i] > Xmax)	x[i] = 2*Xmax - x[i];
			if (y[i] < Ymin)	y[i] = 2*Ymin - y[i];
			if (y[i] > Ymax)	y[i] = 2*Ymax - y[i];
		}
		clock_gettime(CLOCK_MONOTONIC, &time32);
		clock_gettime(CLOCK_MONOTONIC, &time41);
// printf(" joined!\n");
		if (enable) {
			updateWindow();
		}
		// sleep(20);
		clock_gettime(CLOCK_MONOTONIC, &time42);
		printf(" CLOCK_IO_TIME: %lf s\n", (time12.tv_sec - time11.tv_sec) + (time12.tv_nsec - time11.tv_nsec) / 1e9);
		printf(" CLOCK_buildTree_TIME: %lf s\n", (time22.tv_sec - time21.tv_sec) + (time22.tv_nsec - time21.tv_nsec) / 1e9);
		printf(" CLOCK_calculate_TIME: %lf s\n", (time32.tv_sec - time31.tv_sec) + (time32.tv_nsec - time31.tv_nsec) / 1e9);
		printf(" CLOCK_drawWindow_TIME: %lf s\n", (time42.tv_sec - time41.tv_sec) + (time42.tv_nsec - time41.tv_nsec) / 1e9);
	}


    pthread_mutex_destroy(&count_mutex);

    //free
    free(x);
    free(y);
    free(vx);
    free(vy);
    free(Fx);
    free(Fy);
    free(father);
    for (i = 0; i < 5*MAXN; i++)
		free(child[i]);
	free(child);
	free(cnt);
	free(Sm);
	free(Ax);
	free(Ay);
	free(X1);
	free(X2);
	free(Y1);
	free(Y2);

	return 0;
}