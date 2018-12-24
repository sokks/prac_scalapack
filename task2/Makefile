
OPT=	-O2
FLAGS=	-DUSE_DOUBLE

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

start: bin/start
	./bin/start

bin/start: build/start.o lib/smatrix.a Makefile
	mkdir -p bin
	mpicxx -std=c++11 -o $@ $(OPT) build/start.o lib/smatrix.a $(LIBS)

build/start.o: src/start.cpp Makefile
	mkdir -p build
	mpicxx -std=c++11 -o $@ $(OPT) -c src/start.cpp $(FLAGS)


lib/smatrix.a: build/smatrix.o Makefile
	mkdir -p lib
	ar rcs $@ build/smatrix.o

build/smatrix.o: src/smatrix.cpp src/smatrix.h Makefile
	mkdir -p build
	mpicxx -std=c++11 -o $@ $(OPT) -c src/smatrix.cpp $(FLAGS) $(EFLAGS)




clean:
	rm -rf build bin lib doc
