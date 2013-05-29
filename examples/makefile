FC = gfortran

FFLAGS=-O3 -g 
SOURCE = example1.f
TARGET = example1

OBJS1=\
../src/YukawaFMM.o \
../src/rkbesl.o \
../src/ribesl.o \
../src/yukq2d.o


example1: ${OBJS1} example1.o 
	${FC} ${FFLAGS} example1.o -o ${TARGET} \
	${OBJS1} 

ribesl.o:    ../src/ribesl.f
	${FC} ${FFLAGS} -c ../src/ribesl.f

rkbesl.o:    ../src/rkbesl.f
	${FC} ${FFLAGS} -c ../src/rkbesl.f

YukawaFMM.o:  ../src/YukawaFMM.f
	${FC} ${FFLAGS} -c ../src/YukawaFMM.f

example1.o:  example1.f
	${FC} ${FFLAGS} -c example1.f

yukq2d.o:  yukq2d.f
	${FC} ${FFLAGS} -c yukq2d.f


.PHONY: clean
clean:
	rm -rf example1.o
	rm -rf ${OBJS1}
