INCLUDES = -I$(shell root-config --incdir) -I$(NUMIANA_INC)
DK2NU_INCLUDES = -I$(PPFX_DIR)/include -I$(shell root-config --incdir) -I$(BOOSTROOT) -I${DK2NU}/include -I$(NUMIANA_INC)
DEPLIBS=$(shell root-config --libs)

CC	=	g++
COPTS	=	-fPIC -g $(shell root-config --cflags) -std=c++17 -shared -Wno-unused-variable

SRCS = $(shell ls src/*.cxx)
HEADERS = $(shell ls include/*.h | sed -e 's/include\///g')

.PHONY: base flugg dk2nu

all: base flugg dk2nu

base: $(SRCS)
	@echo -e "\033[0;36mBuilding base classes\033[0m"
	if [ ! -d lib ]; then mkdir -p lib; fi
	$(CC) $(COPTS) -shared -o lib/libnumi.so $^ $(DEPLIBS) $(INCLUDES)

flugg: FluggDict.cxx
	@echo -e "\033[0;36mBuilding flugg\033[0m"
	$(CC) $(COPTS) flugg/FluggTree.C -o flugg/FluggTree_C.so $(DEPLIBS) $(INCLUDES)
	$(CC) $(COPTS) dict/FluggDict.cxx flugg/FluggFlux.cc -o lib/FluggFlux_cc.so $(DEPLIBS) $(INCLUDES) -I$(NUMIANA_DIR)/flugg -L$(PWD)/lib -lnumi

FluggDict.cxx:
	@echo -e "\033[0;36mMaking flugg dictionaries\033[0m"
	if [ ! -d dict ]; then mkdir -p dict; fi
	rootcling -f dict/$@ -c $(INCLUDES) -I$(NUMIANA_DIR)/flugg -p $(HEADERS) FluggFlux.h LinkDef.h

dk2nu: Dk2NuDict.cxx
	@echo -e "\033[0;36mBuilding dk2nu\033[0m"
	$(CC) $(COPTS) dict/Dk2NuDict.cxx dk2nu/Dk2NuFlux.cc -o lib/Dk2NuFlux_cc.so $(DEPLIBS) -lEG -L$(PWD)/lib -lnumi -L$(PPFX_DIR)/lib -lppfx -L${DK2NU_LIB} -ldk2nuTree $(DK2NU_INCLUDES) -I$(NUMIANA_DIR)/dk2nu

Dk2NuDict.cxx:
	@echo -e "\033[0;36mMaking dk2nu dictionaries\033[0m"
	if [ ! -d dict ]; then mkdir -p dict; fi
	rootcling -f dict/$@ -c $(DK2NU_INCLUDES) -I$(NUMIANA_DIR)/dk2nu -p $(HEADERS) Dk2NuFlux.h LinkDef.h

clean:
	@echo -e "\033[0;36mCleaning up\033[0m"
	rm -rv lib/*
	rm dict/*Dict.cxx
	rm dict/*.pcm
	rm flugg/*.so
