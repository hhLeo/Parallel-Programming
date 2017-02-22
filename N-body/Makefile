all:Nbody_pthread.out Nbody_baoli.out Nbody_cuda.out Nbody_baoli_cuda.out

# MSC_pthread_static.out:MSC_pthread_static.c
# 	gcc MSC_pthread_static.c -lpthread -lX11 -lrt -o MSC_pthread_static.out

Nbody_pthread.out:Nbody_pthread.c
	g++ Nbody_pthread.c -lpthread -lX11 -lrt -o Nbody_pthread.out

Nbody_baoli.out:Nbody_baoli.c
	g++ Nbody_baoli.c -lpthread -lX11 -lrt -o Nbody_baoli.out

Nbody_cuda.out:Nbody_cuda.cu
	nvcc -arch=sm_20 Nbody_cuda.cu -lX11 -lrt -o Nbody_cuda.out

Nbody_baoli_cuda.out:Nbody_baoli_cuda.cu
	nvcc -arch=sm_20 Nbody_baoli_cuda.cu -lX11 -lrt -o Nbody_baoli_cuda.out

clean:
	rm *.out	
