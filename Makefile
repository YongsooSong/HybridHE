PREFIX=$(HOME)
CC=g++
AR=ar
CFLAGS= -Wall -O3 -std=c++11


LDLIBS = -L/usr/local/lib
LDFLAGS= -I/usr/local/include
LIBS = -lfftw3

INCLUDE= distrib.h params.h FFT.h Scheme.h Database.h Encoding.h

SRC =  distrib.cpp FFT.cpp Scheme.cpp  Database.cpp Encoding.cpp

OBJ= distrib.o FFT.o Scheme.o Database.o Encoding.o

all: libidash2016.a  

#cmd: cmd/gen cmd/enc cmd/nand cmd/dec

install: $(INCLUDE) libidash2016.a
	install $(INCLUDE) $(PREFIX)/include
	install libidash2016.a $(PREFIX)/lib
new:
	make clean
	make all
	g++ -I/usr/local/include main.cpp libidash2016.a -o foo -O3 -L/usr/local/lib -lfftw3
	
clean:
	rm *.o libidash2016.a  || echo nothing to clean

query1:
	./foo 1 164129713 C T RCV000015246_10000.txt

query2:
	./foo 1 161276680 A T RCV000015246_100000.txt

query3:
	./foo 1 164129713 A T RCV000015246_100000.txt

test:
	./foo 1 161150579 G A RCV000015246_test.txt
	
obj: $(OBJ)	

libidash2016.a: $(OBJ)
	$(AR) -q libidash2016.a  $(OBJ)

%.o: %.cpp $(INCLUDE)
	 $(CC) $(CFLAGS) -c $< $(LDFLAGS)
