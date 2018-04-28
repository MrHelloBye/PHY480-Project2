jacobi: main.o eigen.o
	g++ main.o eigen.o -o jacobi

main.o: main.cpp
	g++ main.cpp -c -std=c++11 -Ofast

eigen.o: eigen.cpp eigen.h
	g++ eigen.cpp -c -std=c++11 -Ofast

unit_tests: main_test.cpp test_functions.h eigen.cpp eigen.h
	g++ eigen.cpp -c -std=c++11 -Ofast
	g++ main_test.cpp -c -std=c++11 -Ofast
	g++ main_test.o eigen.o -o unit_tests
	./unit_tests

clean:
	rm *.o
