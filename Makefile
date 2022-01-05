build:
	mpicc main.c libbmp.c leader.c common.c follower.c -lgmp -lm -g 

run:
	mpirun --mca opal_warn_on_missing_libcuda 0 --oversubscribe -n 5 a.out -i 1000 -z 3 -s 256 -e 

debug:
	mpirun --mca opal_warn_on_missing_libcuda 0 --oversubscribe -n 1 gdb a.out : -n 4 a.out -i 50 -z 2 -s 256

clean:
	rm -f *.bmp *.bin a.out
