FC = gfortran

FFLAGS=-O3 -g -s

OBJS1=\
../src/YukawaFMM.o \
../src/rkbesl.o \
../src/ribesl.o \
../src/yukq2d.o

OBJS2=\
../src/btos_col.o \
../src/btos_double.o \
../src/btos_single.o \
../src/col_double.o \
../src/col_loc.o \
../src/col_single.o \
../src/forced.o \
../src/YukawaFMM.o \
../src/rkbesl.o \
../src/ribesl.o \
../src/stob_col.o \
../src/stob_double.o \
../src/stob_single.o \
../src/tab_cheby.o \
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

tab_cheby.o: tab_cheby.f
	${FC} ${FFLAGS} -c tab_cheby.f

btos_col.o: btos_col.f
	${FC} ${FFLAGS} -c btos_col.f

btos_double.o: btos_double.f
	${FC} ${FFLAGS} -c btos_double.f

btos_single.o: btos_single.f
	${FC} ${FFLAGS} -c btos_single.f

col_double.o: col_double.f
	${FC} ${FFLAGS} -c col_double.f

col_loc.o: col_loc.f
	${FC} ${FFLAGS} -c col_loc.f

col_single.o: col_single.f
	${FC} ${FFLAGS} -c col_single.f

stob_col.o: stob_col.f
	${FC} ${FFLAGS} -c stob_col.f

stob_double.o: stob_double.f
	${FC} ${FFLAGS} -c stob_double.f

stob_single.o: stob_single.f
	${FC} ${FFLAGS} -c stob_single.f


.PHONY: clean
clean:
	rm -rf example1.o
	rm -rf example2.o
	rm -rf ${OBJS1}
