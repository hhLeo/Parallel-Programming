all:kn.out kn_fstream.out

kn.out:kn.cu
	nvcc -arch=sm_20 kn.cu -o kn.out

kn_fstream.out:kn_fstream.cu
	nvcc -arch=sm_20 kn_fstream.cu -o kn_fstream.out

clean:
	rm *.out