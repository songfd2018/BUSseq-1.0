CC=g++

all: BUSseq BUSseq_inference

# link
BUSseq: BUSseq.o rngstream.o
	$(CC) -fopenmp BUSseq.o rngstream.o -o BUSseq

BUSseq_inference: BUSseq_inference.o rngstream.o
	$(CC) -fopenmp BUSseq_inference.o rngstream.o -o BUSseq_inference

# compile
BUSseq.o: BUSseq.cpp
	$(CC) -c -fopenmp BUSseq.cpp -w -std=c++11

BUSseq_inference.o: BUSseq_inference.cpp
	$(CC) -c -fopenmp BUSseq_inference.cpp -w -std=c++11

rngstream.o: rngstream.cpp
	$(CC) -c -fopenmp -lgomp rngstream.cpp

# clean all files
clean:
	rm -f *.o;