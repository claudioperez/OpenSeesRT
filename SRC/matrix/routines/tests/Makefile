cmx = ../inv2.o ../inv3.o ../inv4.o ../inv5.o ../inv6.o
arch = #-march=haswell
test: test_inverse.c $(cmx) Makefile
	cc test_inverse.c inv[2-6].o -llapack -lblas -lm -o test

%.o: %.c Makefile cmx.h
	cc -Ofast -c $< -llapack -lblas -lm -o $@ \
	 -fno-math-errno -fno-signaling-nans -fno-trapping-math \
	 -fassociative-math -ffast-math -mtune=native -mfpmath=sse $(arch)
