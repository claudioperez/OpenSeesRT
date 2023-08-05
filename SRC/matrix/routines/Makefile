
test: test_inverse.c inv6.o Makefile
	cc test_inverse.c inv6.o -llapack -lblas -lm -o test

%.o: %.c Makefile
	cc -Ofast -c $< -llapack -lblas -lm -o $@ \
	 -fno-math-errno -fno-signaling-nans -fno-trapping-math \
	 -fassociative-math -ffast-math -march=haswell -mtune=native -mfpmath=sse
