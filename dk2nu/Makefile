INCLUDES = -I$(PPFX_DIR)/include -I$(shell root-config --incdir) -I$(BOOSTROOT) -I${DK2NU}/include
DEPLIBS=$(shell root-config --libs) -lEG

CC	=	g++
COPTS	=	-fPIC -g $(shell root-config --cflags) -std=c++17 -shared

dk2nu:
	$(CC) $(COPTS) Dk2NuFlux.cc -o Dk2NuFlux_cc.so $(DEPLIBS) -L$(PPFX_DIR)/lib -lppfx -L${DK2NU_LIB} -ldk2nuTree $(INCLUDES)
