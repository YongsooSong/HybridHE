Our code requires the FFTW 3 library for polynomial arithmetic which is available at http://www.fftw.org/download.html, and a c++ compiler.

$make clean
$make all

This will build the code. Then the following should work:

$g++ main.cpp libidash2016.a -o foo -O3 -lfftw3

for compiling the code and making your program foo.c
The program will run with a query consisting of (#Chr, Position, REF, ALT, VCF filename).
For example, you can write:

$./foo 1,1,1 161276680,161276217,161276672 A,G,C T,C,T RCV000015246_100000.txt

You can write ‘E’ for insertion or deletion query.
(e.g. $./foo 1 73934717 E T RCV000015246_1000000.txt)

Our program compares the first 20 characters of REF and ALT (of query and database) for efficiency, 
and outputs the matching result.

“nPoly” is the number of RLWE ciphertexts for database.
FFT setup prepares polynomial arithmetic (independent of user database or query), 
so the server may precompute it.