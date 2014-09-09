# Apparently Apple's compiler stores string constants in char[1] buffers.
# Boost treats char[1]'s as empty paths and disregards them assuming that a char[1] is only large
# enough to hold the terminal null character.
# Because of this, we cannot use Apple's LVM compiler when working with boost::filesystem::path and is_directory(SOME_PATH)
# Below we default to the custom build G++ version in /usr/local/bin

# Select compiler
COMPILER = GNU


OS := $(shell uname)

ifeq ($(COMPILER),GNU)

# Find out what OS we're working on.
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  # Custom built GNU compiler on OS X since what ships with Xcode currently lacks
  # OMP support
  CXX	   = /usr/local/bin/g++
else
  CXX	   = /usr/bin/g++
endif

endif



# Makefile is giving an error.  This is a temporary fix.
#g++ *.cpp -lboost_thread -lboost_system -lboost_filesystem
# CFLAGS for running
CFLAGS = -Wall -mtune=native  -m64 -msse4.2 -O3

# CFLAGS for debugging
#CFLAGS = -Wall -O0 -g

RM = rm -rf
LIBS =-lboost_thread -lboost_filesystem -lboost_system

SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)



.cpp.o:
	 $(CXX) -c $(CFLAGS) $*.cpp


all : speckle-boost 


speckle-boost: $(OBJS)
	 $(CXX)  -o  $@ $(OBJS) $(CFLAGS) $(LIBS)

depend:
	mkdep ${CFLAGS} ${SRCS}

-include .depend

clean::
	 $(RM) speckle-boost*
	 $(RM) *.o
	 $(RM) *.depend

