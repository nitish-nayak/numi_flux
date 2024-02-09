INCLUDES = -I$(shell root-config --incdir) -I$(NUMIANA_INC)
DK2NU_INCLUDES = -I$(PPFX_DIR)/include -I$(shell root-config --incdir) -I$(BOOSTROOT) -I${DK2NU}/include -I$(NUMIANA_INC)
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
	$(CC) $(COPTS) dict/FluggDict.cxx flugg/FluggFlux.cc -o lib/FluggFlux_cc.so $(DEPLIBS) $(INCLUDES) -I$(NUMIANA_DIR)/flugg -L$(PWD)/lib -lnumi

FluggDict.cxx:
	if [ ! -d dict ]; then mkdir -p dict; fi
	rootcling -f dict/$@ -c $(INCLUDES) -I$(NUMIANA_DIR)/flugg -p $(HEADERS) FluggFlux.h LinkDef.h

dk2nu: Dk2NuDict.cxx
	$(CC) $(COPTS) dict/Dk2NuDict.cxx dk2nu/Dk2NuFlux.cc -o lib/Dk2NuFlux_cc.so $(DEPLIBS) -lEG -L$(PWD)/lib -lnumi -L$(PPFX_DIR)/lib -lppfx -L${DK2NU_LIB} -ldk2nuTree $(DK2NU_INCLUDES) -I$(NUMIANA_DIR)/dk2nu

Dk2NuDict.cxx:
	if [ ! -d dict ]; then mkdir -p dict; fi
	rootcling -f dict/$@ -c $(DK2NU_INCLUDES) -I$(NUMIANA_DIR)/dk2nu -p $(HEADERS) Dk2NuFlux.h LinkDef.h

clean:
	rm -rv lib/*
	rm dict/*Dict.cxx
	rm dict/*.pcm
	rm flugg/*.so
