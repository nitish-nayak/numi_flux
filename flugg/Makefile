INCLUDES = -I$(PWD)
DEPLIBS=$(shell root-config --libs)

CC	=	g++
COPTS	=	-fPIC -g $(shell root-config --cflags) -std=c++17 -shared

flugg:
	$(CC) $(COPTS) FluggTree.C -o FluggTree_C.so $(DEPLIBS) $(INCLUDES)
	$(CC) $(COPTS) FluggFlux.cc -o FluggFlux_cc.so $(DEPLIBS) $(INCLUDES)

clean:
	rm *.so
