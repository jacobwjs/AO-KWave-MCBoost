# Makefile is giving an error.  This is a temporary fix.
#g++ *.cpp -lboost_thread -lboost_system -lboost_filesystem
# CFLAGS for running
CFLAGS = -Wall -mtune=native  -m64 -msse4.2 -O3

# CFLAGS for debugging
#CFLAGS = -Wall -O0 -g

CC = g++
RM = rm -rf
LIBS =-lboost_thread -lboost_filesystem -lboost_system

SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)



.cpp.o:
	 $(CC) -c $(CFLAGS) $*.cpp 


all : speckle-boost 


speckle-boost: $(OBJS)
	 $(CC)  -o  $@ $(OBJS) $(CFLAGS) $(LIBS)

depend:
	mkdep ${CFLAGS} ${SRCS}

-include .depend

clean::
	 $(RM) speckle-boost*
	 $(RM) *.o
	 $(RM) *.depend

