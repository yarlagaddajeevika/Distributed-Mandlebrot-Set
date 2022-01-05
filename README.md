# 474-Project-2
Project 2: Implement a distributed algorithm using MPI

## Distributed Mandelbrot Set

This project calculates the mandelbrot set with 
arbitrary percision floating point numbers and distributes 
the work using MPI. 

### Group members

Dalton Caron dpcaron@csu.fullerton.edu
Jeevika Yarlagadda jeevikayarlagadda@csu.fullerton.edu

### Compilation

The user must have GCC, libgmp, and make installed to compile the program. The user's compiler  
must comply to POSIX standards. To compile the program, use the Makefile as shown below. 

```bash
make build
```

This command produces the a.out executable file. 

### How to use the program

This program produces bitmap image files and binary files storing the parameters used to produce 
the images. There are a variety of command line options:

- -e: save the mandelbrot bmp every zoom iteration of the algorithm
- -f [file]: load parameters from a binary save file 
- -z [number]: the amount of times the algorithm will recompute the set with a deeper zoom
- -s [number]: set the size of the output image. Must be a power of 2. 
- -i [number]: the amount of times to iterate over the formula f(z). Must be larger the deeper the zoom. 

An example command is shown below. 

```bash
mpirun --mca opal_warn_on_missing_libcuda 0 --oversubscribe -n 5 a.out -i 1000 -z 3 -s 256  
```

Depending on the parameters, the algorithm may execute in a few seconds or few hours. The default parameters in the 
Makefile run in a few seconds on modern machines. The Makefile default is executed using the following command.

```bash
make run
```

There are debugging presets, but please disregard them unless you are a developer. 

### Pseudocode

The pseudocode for this project is located in `report/report.pdf` in the root of this repository. 
