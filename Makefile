#! /usr/bin/make

CC=mpic++
CXXFLAGS=-std=c++0x

# you need to change this
DASHROOT=${HOME}/work/dash

	INC=-I${DASHROOT}/dart-impl/mpi/dart-mpi/ -I${DASHROOT}/dart-if/v3 -I${DASHROOT}/dash/dash-lib/
LIBS=${DASHROOT}/dash/dash-lib/libdash.a ${DASHROOT}/dart-impl/mpi/dart-mpi/libdart.a

SOURCES= heat2d.cpp

OBJS=$(SOURCES:%.cpp=%.o)

.PHONY: clean
.PHONY: all

%.o:%.cpp
	$(CC) $(CXXFLAGS) $(INC) -o $@ $< $(LIBS)

all: $(OBJS)

clean:
	rm -rf *.o
