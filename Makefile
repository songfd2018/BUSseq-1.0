CC=g++

all: BUSseq_od BUSseq_inference_od

# link
BUSseq_od: BUSseq_od.o rngstream.o
	$(CC) -fopenmp BUSseq_od.o rngstream.o -o BUSseq_od

BUSseq_inference_od: BUSseq_inference_od.o rngstream.o
	$(CC) -fopenmp BUSseq_inference_od.o rngstream.o -o BUSseq_inference_od

# compile
BUSseq_od.o: BUSseq_od.cpp
	$(CC) -c -fopenmp BUSseq_od.cpp -w -std=c++11

BUSseq_inference_od.o: BUSseq_inference_od.cpp
	$(CC) -c -fopenmp BUSseq_inference_od.cpp -w -std=c++11

rngstream.o: rngstream.cpp
	$(CC) -c -fopenmp -lgomp rngstream.cpp

# clean all files
clean:
	rm -f *.o;