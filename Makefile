INCLUDES = -I./include -I$(shell root-config --incdir)
DEPLIBS=$(shell root-config --libs)

CC	=	g++
COPTS	=	-fPIC -g $(shell root-config --cflags) -std=c++17 -shared

SRCS = $(shell ls src/*.cxx)

base: $(SRCS)
	if [ ! -d lib ]; then mkdir -p lib; fi

	$(CC) $(COPTS) -shared -o lib/libnumi.so $^ $(DEPLIBS) $(INCLUDES)

clean:
	rm -rv lib/*
