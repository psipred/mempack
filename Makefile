# Makefile for MEMPACK-SVM and SVM_CLASSIFY

CC=gcc
CPP=g++
LD=gcc
CFLAGS=-O3
LFLAGS=-O3
LIBS=-lm
INC=/usr/include/
BOOST=/opt/boost_1_37_0
MKDIR=mkdir

all: create_input create_bin create_output svm_classify kk_plot

clean:
	rm -f bin/svm_classify
	rm -f bin/kk_plot
	rm -f src/svm_classify.o
	rm -f src/svm_common.o

create_input:
	$(MKDIR) -p input/

create_bin:
	$(MKDIR) -p bin/

create_output:
	$(MKDIR) -p output/

svm_common.o: src/svm_common.c src/svm_common.h src/kernel.h
	$(GCC) -c $(CFLAGS) src/svm_common.c -o src/svm_common.o

svm_classify.o: svm_classify.c svm_common.h kernel.h
	$(GCC) -c $(CFLAGS) src/svm_classify.c -o src/svm_classify.o

svm_classify: src/svm_classify.o src/svm_common.o
	$(LD) $(LFLAGS) src/svm_classify.o src/svm_common.o -o bin/svm_classify $(LIBS)

kk_plot: src/draw_graphs.cpp src/globals.cpp src/paramopt.c
	$(CPP) --std=c++11 -Wno-write-strings -Wno-deprecated -I$(BOOST) -I$(INC) $(LIBS) -O2 src/draw_graphs.cpp src/globals.cpp src/paramopt.c -o bin/kk_plot
