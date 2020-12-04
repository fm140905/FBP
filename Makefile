# Makefile for main

# *****************************************************
# Variables to control Makefile operation

CXX = g++
CXXFLAGS = -Wall -g -std=c++11
# CXXFLAGS = -static -static-libgcc -static-libstdc++

VPATH = source

OUT = bin/


OBJECTS = $(OUT)main.o $(OUT)cpimaging.o $(OUT)getInput.o



LIBS1 = `root-config --cflags --glibs`
LIBS2 = `root-config --cflags --glibs`

# ****************************************************
# Targets needed to bring the executable up to date

$(OUT)main: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LIBS1)

	# The main target

$(OUT)main.o: main.cc getInput.h cpimaging.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBS1)

	# The main.o target

$(OUT)getInput.o: getInput.cc getInput.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBS2)

	# The getInput.o target

$(OUT)cpimaging.o: cpimaging.cc cpimaging.h getInput.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBS2)

	# The cpimaging.o target

.PHONY: clean

clean :
	rm $(OUT)*.o