

classify : classify.cpp gzstream/gzstream.C gzstream/gzstream.h kmer/kmer.h
	g++ -g -c  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11 classify.cpp gzstream.o -lz -lpthread -o classify
