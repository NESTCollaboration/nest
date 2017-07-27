
all: fastNEST

fastNEST: libNEST.cpp LUX.hh
	g++ -Ofast -std=c++0x libNEST.cpp -o fastNEST.exe

clean:
	rm -f fastNEST*
