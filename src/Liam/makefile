ho: main_ho.o eigen.o
	g++ main_ho.o eigen.o -o ho

qd: main_qd.o eigen.o
	g++ main_qd.o eigen.o -o qd

main_ho.o: main_ho.cpp
	g++ main_ho.cpp -c -std=c++11 -Ofast

main_qd.o: main_qd.cpp
	g++ main_qd.cpp -c -std=c++11 -Ofast

eigen.o: eigen.cpp eigen.h
	g++ eigen.cpp -c -std=c++11 -Ofast

unit_tests: main_test.cpp test_functions.h eigen.cpp eigen.h
	g++ eigen.cpp -g -c -std=c++11 -Ofast
	g++ main_test.cpp -g -c -std=c++11 -Ofast
	g++ main_test.o eigen.o -o unit_tests

clean:
	rm *.o
