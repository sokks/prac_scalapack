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


all: bin

bin: bin/test

lib: lib/smatrix.a lib/soutput.a

build/test.o: src/test.cpp Makefile
	mkdir -p build
	mpicxx -o $@ $(OPT) -c src/test.cpp $(FLAGS)

build/smatrix.o: src/smatrix.cpp Makefile
	mkdir -p build
	mpicxx -o $@ $(OPT) -c src/smatrix.cpp $(FLAGS) $(EFLAGS)

lib/smatrix.a: build/smatrix.o Makefile
	mkdir -p lib
	ar rcs $@ build/smatrix.o
bin/test: build/test.o lib/smatrix.a Makefile
	mkdir -p bin
	mpicxx -o $@ $(OPT) build/test.o lib/smatrix.a $(LIBS)

test: bin/test
	mpiexec -np 4 ./bin/test

clean:
	rm -rf build
	find . -name \*~ -delete

clear:
	rm -rf bin lib doc

doc: src Doxyfile Makefile
	rm -rf doc
	doxygen Doxyfile

doxywizard:
	doxywizard Doxyfile

archive: clean
	git archive --format tar master | xz -9 > matrix_library.tar.xz
