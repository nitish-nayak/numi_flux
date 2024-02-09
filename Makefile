INCLUDES = -I./include -I$(shell root-config --incdir)
DK2NU_INCLUDES = -I./include -I$(PPFX_DIR)/include -I$(shell root-config --incdir) -I$(BOOSTROOT) -I${DK2NU}/include
DEPLIBS=$(shell root-config --libs)

CC	=	g++
COPTS	=	-fPIC -g $(shell root-config --cflags) -std=c++17 -shared

SRCS = $(shell ls src/*.cxx)
HEADERS = $(shell ls include/*.h | sed -e 's/include\///g')

.PHONY: base flugg dk2nu

all: base flugg dk2nu

base: $(SRCS)
	if [ ! -d lib ]; then mkdir -p lib; fi
	$(CC) $(COPTS) -shared -o lib/libnumi.so $^ $(DEPLIBS) $(INCLUDES)

flugg: FluggDict.cxx
	$(CC) $(COPTS) flugg/FluggTree.C -o flugg/FluggTree_C.so $(DEPLIBS) $(INCLUDES)
	$(CC) $(COPTS) FluggDict.cxx flugg/FluggFlux.cc -o FluggFlux_cc.so $(DEPLIBS) $(INCLUDES) -I./flugg -L$(PWD)/lib -lnumi

FluggDict.cxx:
	rootcling -f $@ -c $(INCLUDES) -I./flugg -p $(HEADERS) FluggFlux.h LinkDef.h

dk2nu: Dk2NuDict.cxx
	$(CC) $(COPTS) Dk2NuDict.cxx dk2nu/Dk2NuFlux.cc -o Dk2NuFlux_cc.so $(DEPLIBS) -lEG -L$(PWD)/lib -lnumi -L$(PPFX_DIR)/lib -lppfx -L${DK2NU_LIB} -ldk2nuTree $(DK2NU_INCLUDES) -I./dk2nu

Dk2NuDict.cxx:
	rootcling -f $@ -c $(DK2NU_INCLUDES) -I./dk2nu -p $(HEADERS) Dk2NuFlux.h LinkDef.h

clean:
	rm -rv lib/*
	rm *Dict.cxx
	rm *.pcm
	rm *.so
	rm flugg/*.so
	rm dk2nu/*.so
