INCLUDES = -I./include -I$(shell root-config --incdir)
DEPLIBS=$(shell root-config --libs)

CC	=	g++
COPTS	=	-fPIC -g $(shell root-config --cflags) -std=c++17 -shared

SRCS = $(shell ls src/*.cxx)
HEADERS = $(shell ls include/*.h | sed -e 's/include\///g')

.PHONY: base flugg dk2nu

base: $(SRCS)
	if [ ! -d lib ]; then mkdir -p lib; fi
	$(CC) $(COPTS) -shared -o lib/libnumi.so $^ $(DEPLIBS) $(INCLUDES)

flugg: FluggDict.cxx
	$(CC) $(COPTS) flugg/FluggTree.C -o flugg/FluggTree_C.so $(DEPLIBS) $(INCLUDES)
	$(CC) $(COPTS) FluggDict.cxx flugg/FluggFlux.cc -o FluggFlux_cc.so $(DEPLIBS) $(INCLUDES) -I./flugg -L$(PWD)/lib -lnumi

FluggDict.cxx:
	rootcling -f $@ -c $(INCLUDES) -I./flugg -p $(HEADERS) FluggFlux.h LinkDef.h

clean:
	rm -rv lib/*
	rm *Dict.cxx
	rm *.pcm
	rm *.so
	rm flugg/*.so
