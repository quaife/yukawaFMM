FC = gfortran

FFLAGS=-O3 -g 

OBJS1=\
../src/YukawaFMM.o \
../src/rkbesl.o \
../src/ribesl.o \
../src/yukq2d.o

OBJS2=\
../src/forced.o \
../src/YukawaFMM.o \
../src/rkbesl.o \
../src/ribesl.o \
../src/yukq2d.o

example1: ${OBJS1} example1.o 
	${FC} ${FFLAGS} example1.o -o example1 \
	${OBJS1} 

example2: ${OBJS2} example2.o 
	${FC} ${FFLAGS} example2.o -o example2 \
	${OBJS2} 

ribesl.o:    ../src/ribesl.f
	${FC} ${FFLAGS} -c ../src/ribesl.f

rkbesl.o:    ../src/rkbesl.f
	${FC} ${FFLAGS} -c ../src/rkbesl.f

YukawaFMM.o:  ../src/YukawaFMM.f
	${FC} ${FFLAGS} -c ../src/YukawaFMM.f

example1.o:  example1.f
	${FC} ${FFLAGS} -c example1.f

example2.o:  example2.f
	${FC} ${FFLAGS} -c example2.f

yukq2d.o:  yukq2d.f
	${FC} ${FFLAGS} -c yukq2d.f

forced.o:  forced.f
	${FC} ${FFLAGS} -c forced.f


.PHONY: clean
clean:
	rm -rf example1.o
	rm -rf example2.o
	rm -rf ${OBJS1}
