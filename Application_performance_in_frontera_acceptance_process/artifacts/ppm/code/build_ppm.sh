#build copen/close/  program
mpicc -c -g -mcmodel=medium cio.c

#build main program
mpifort -mcmodel=medium -qopenmp -xCORE-AVX512 -O3 -o ppm2f_5120_33792_tasks_4threads_ompi PPM2F-tp3-5-1-12-ICF-short-loops.F cio.o 
