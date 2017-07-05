
all: fastNEST

fastNEST: libNEST.cpp LUX.hh
	g++ -Ofast libNEST.cpp -o fastNEST.exe

clean:
	rm -f fastNEST*
