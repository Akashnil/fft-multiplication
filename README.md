# fft-multiplication
A small python implementation for large number multiplication using fft. We test and demonstrate
~ O(n log n) complexity. It is not Schönhage–Strassen_algorithm in its full complexity. It is
also not optimized for performance, kept simple enough to understand.

The simpler implementation uses complex numbers: in complex-arithmetic.py

This implementation is faster for base 2 multiplication, but for larger bases (or bery large number
of digits) gives wrong answer because of finite precision floating point calculations

The more nuanced implementation uses a finite field: in galois-field-arithmetic.py

This implementation can be extended to larger bases and number of digits. For binary numbers,
better performance is achieved by grouping digits and working in base 2^k for small k, eg. 5.
