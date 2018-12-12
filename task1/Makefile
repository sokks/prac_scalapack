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

test_mwriter: bin/mwriter_test
	./bin/mwriter_test data/testm.txt data/testdm.txt data/testm.dat data/testdm.dat

test_smatrix: bin/test
	mpiexec -np 1 ./bin/test 3 ./data/R 0.5 ./data/H 10

lib/mwriter.a: build/mwriter.o
	mkdir -p lib
	ar rcs $@ build/mwriter.o

lib/smatrix.a: build/smatrix.o Makefile
	mkdir -p lib
	ar rcs $@ build/smatrix.o

build/mwriter.o: src/mwriter/mwriter.cpp Makefile
	mkdir -p build
	g++ -std=c++11 -o $@ -c src/mwriter/mwriter.cpp

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

gen_diag: bin/gen_diag
	# последний параметр -- seed (если не задан берется текущее время) 
	./bin/gen_diag data/ro_diag_100.dat 100 100

bin/gen_diag: lib/mwriter.a src/gen_ro_diag.cpp
	mkdir -p bin
	g++ -std=c++11 -o $@ src/gen_ro_diag.cpp lib/mwriter.a


clean:
	rm -rf build bin lib doc

