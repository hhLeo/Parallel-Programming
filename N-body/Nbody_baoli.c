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

//Window
Display *display;
Window window;
int screen;
GC gc;
XGCValues values;
int width, height;	//window (= Length)

void drawWindow() {
	printf("drawWindow init\n");
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
	printf("update...\n");
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

int main(int argc, char const *argv[])
{
	//malloc
	int i, j;
	x = (double*)malloc(sizeof(double)*MAXN);
	y = (double*)malloc(sizeof(double)*MAXN);
	vx = (double*)malloc(sizeof(double)*MAXN);
	vy = (double*)malloc(sizeof(double)*MAXN);
	Fx = (double*)malloc(sizeof(double)*MAXN);
	Fy = (double*)malloc(sizeof(double)*MAXN);

	printf(" 1. malloc over\n");

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
	Xmin -= 0.2 * dx;
	Xmax += 0.2 * dx;
	Ymin -= 0.2 * dy;
	Ymax += 0.2 * dy;

	if (strcmp(argv[7], "enable") == 0)	enable = 1;
	else	enable = 0;
	printf("enable = %d\n", enable);
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
	printf(" 2. ReadIn over\n");

	printf(" 3. while T start...\n");
	//main
	while(T--) {
		printf(" *** T = %d ***\n", T);

//baoli
		for(i = 0; i < N; i++) {
			Fx[i] = Fy[i] = 0;
			for(j = 0; j < N; j++) {
				if(i == j)continue;
				double r2 = (x[j]-x[i])*(x[j]-x[i])+(y[j]-y[i])*(y[j]-y[i]);
				double r3 = r2*sqrt(r2);
				if(sqrt(r2) < 1e-10)	continue;
				Fx[i] += G*m*m*(x[j]-x[i])/r3;
				Fy[i] += G*m*m*(y[j]-y[i])/r3;
			}
		}
		for(i = 0; i < N; i++) {
			vx[i] += Fx[i]*t/m;
			vy[i] += Fy[i]*t/m;
			x[i] += vx[i]*t;
			y[i] += vy[i]*t;
				// if (x[i] < Xmin)	x[i] = 2*Xmin - x[i];
				// if (x[i] > Xmax)	x[i] = 2*Xmax - x[i];
				// if (y[i] < Ymin)	y[i] = 2*Ymin - y[i];
				// if (y[i] > Ymax)	y[i] = 2*Ymax - y[i];
			
		}
	
		if (enable) {
			updateWindow();
		}
		// sleep(100);
	}

	//time
/*
	struct timespec time1 = {0, 0};
	struct timespec time2 = {0, 0};
	clock_gettime(CLOCK_MONOTONIC, &time1);
	//do sth...
	clock_gettime(CLOCK_MONOTONIC, &time2);
	printf("CLOCK_MONOTONIC_DELTA_TIME: %lf ms\n", (time2.tv_sec - time1.tv_sec) * 1000 + (time2.tv_nsec - time1.tv_nsec) / 1e6);
*/

	return 0;
}