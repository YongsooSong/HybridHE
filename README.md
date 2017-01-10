# HybridHE

Our code requires the FFTW 3 library for polynomial arithmetic which is available at http://www.fftw.org/download.html, and a c++ compiler.

$make clean
$make all

This will build the code. Then the following should work:

$g++ main.cpp libidash2016.a -o foo -O3 -lfftw3

for compiling the code and making your program foo.c
The program will run with a query consisting of (#Chr, Position, REF, ALT, VCF filename).
For example, you can write:

$./foo 1 161276680 A T RCV000015246_100000.txt

You can write ‘E’ for insertion or deletion query.

It will compare the first several characters of SNPs for efficiency, and outputs the matching result.
