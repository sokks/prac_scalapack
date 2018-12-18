OPT=	-O2
FLAGS=	-DUSE_COMPLEX -DUSE_DOUBLE

ifneq (,$(findstring lomonosov,$(shell hostname)))
	LIBS=	-lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blacs_intelmpi_lp64 \
		-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
	EFLAGS=
else
	LIBPATH=`pwd`/../scalapack-2.0.2/lib/
	LIBPATH= /usr/local/lib/
	LIBS=	$(LIBPATH)libscalapack.a \
		-lpthread -lblas -llapack -lm -lgfortran
	EFLAGS=
endif



test: test_smatrix



test_smatrix: bin/test
	mpiexec -np 1 ./bin/test 3 ./data/R 0.5 ./data/H 10


lib/smatrix.a: build/smatrix.o Makefile
	mkdir -p lib
	ar rcs $@ build/smatrix.o


build/smatrix.o: src/matrix/smatrix.cpp src/matrix/smatrix.h src/matrix/scalapack.h Makefile
	mkdir -p build
	mpicxx -o $@ $(OPT) -c src/matrix/smatrix.cpp $(FLAGS) $(EFLAGS)

build/test.o: src/test.cpp Makefile
	mkdir -p build
	mpicxx -o $@ $(OPT) -c src/test.cpp $(FLAGS)


bin/mwriter_test: lib/mwriter.a src/mwriter/mwriter_test.cpp
	mkdir -p bin
	g++ -std=c++11 -o $@ src/mwriter/mwriter_test.cpp lib/mwriter.a

bin/test: build/test.o lib/smatrix.a Makefile
	mkdir -p bin
	mpicxx -o $@ $(OPT) build/test.o lib/smatrix.a $(LIBS)





clean:
	rm -rf build bin lib doc



test_mwriter: bin/mwriter_test
	./bin/mwriter_test data/testm.txt data/testdm.txt data/testm.dat data/testdm.dat

lib/mwriter.a: build/mwriter.o
	mkdir -p lib
	ar rcs $@ build/mwriter.o

build/mwriter.o: src/mwriter/mwriter.cpp Makefile
	mkdir -p build
	g++ -std=c++11 -o $@ -c src/mwriter/mwriter.cpp

bin/gen_H: src/gen_H.cpp lib/mwriter.a Makefile
	mkdir -p bin
	g++ -std=c++11 -o $@ src/gen_H.cpp lib/mwriter.a

gen_H: bin/gen_H
	./bin/gen_H data/H1 8 100

bin/gen_ro_full: src/gen_ro_full.cpp lib/mwriter.a Makefile
	mkdir -p bin
	g++ -std=c++11 -o $@ src/gen_ro_full.cpp lib/mwriter.a

gen_ro_full: bin/gen_ro_full
	./bin/gen_ro_full data/R1 8 100

bin/gen_ro_diag: src/gen_ro_diag.cpp lib/mwriter.a Makefile
	mkdir -p bin
	g++ -std=c++11 -o $@ src/gen_ro_diag.cpp lib/mwriter.a

gen_ro_diag: bin/gen_ro_diag
	./bin/gen_ro_diag data/R2 8 0

gen_diag: bin/gen_diag
	# последний параметр -- seed (если не задан берется текущее время) 
	./bin/gen_diag data/ro_diag_100.dat 100 0

bin/gen_diag: lib/mwriter.a src/gen_ro_diag.cpp Makefile
	mkdir -p bin
	g++ -std=c++11 -o $@ src/gen_ro_diag.cpp lib/mwriter.a
