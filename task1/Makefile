test: test_mwriter


test_mwriter: bin/mwriter_test
	./bin/mwriter_test 5 data/testm.in

lib/mwriter.a: build/mwriter.o
	mkdir -p lib
	ar rcs $@ build/mwriter.o

build/mwriter.o: src/mwriter/mwriter.cpp Makefile
	mkdir -p build
	g++ -o $@ -c src/mwriter/mwriter.cpp

bin/mwriter_test: lib/mwriter.a src/mwriter/mwriter_test.cpp
	mkdir -p bin
	g++ -o $@ src/mwriter/mwriter_test.cpp lib/mwriter.a

clean:
	rm -rf build bin lib doc

