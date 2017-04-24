# We use a finite field Z[N] with a suitable root of unity to carry out the FFT

# forward fft using cooley tukey algorithm
def fft(sequence, root, modulo):
    if (len(sequence) == 1):
        return
    even = sequence[::2]
    odd = sequence[1::2]
    r = root * root
    r %= modulo
    fft(even, r, modulo)
    fft(odd, r, modulo)
    r = 1
    m = len(sequence)
    for i in range(m):
        j = i % (m/2)
        sequence[i] = even[j] + odd[j] * r
        r *= root
        sequence[i] %= modulo
        r %= modulo

# inverse fft
def ifft(sequence, root, inverse, modulo):
    m = len(sequence)
    for i in range(1, (m+1)/2):
        sequence[m-i], sequence[i] = sequence[i], sequence[m-i]
    fft(sequence, root, modulo)
    for i in range(m):
        sequence[i] = (sequence[i] * inverse) % modulo

# returns N, r where r has order 2^k modulo N
# N is likely to be prime, but not neccessary.
def find_ring(k):
    # We are looking for N = 2^k * t + 1 and r such that
    # r has order 2^k modulo N, i.e. r^2^(k-1) = -1
    t = 1
    while True:
        N = (t<<k) + 1
        for r in range(2, t+2):
            s = r
            for i in range(k-1):
                s *= s
                s %= N
            if N-s == 1:
                return N, r
        t += 1

# max is the maximum signal
def convolve(x, y, max = 1):
    n = len(x)
    k = len(bin(n))-3
    w = len(bin(max * max))-2
    # find a root or order 2^(k+w)
    modulo, root = find_ring(k+w)
    # root^2^w has order 2^k, which is what we need
    for i in range(w):
        root *= root
        root %= modulo
    inverse = modulo - int((modulo-1)/n)
    fft(x, root, modulo)
    fft(y, root, modulo)
    z = [x[i] * y[i] for i in range(len(x))]
    ifft(z, root, inverse, modulo)
    return z

# x, y are digit arrays (low to high order) in base b
# X = x[0] + x[1]b + x[1]b^2 + ...
def multiply_digits(x, y, b):
    # round up digits to power of 2 with the later half all zeros
    n = 1 << (len(bin(max(len(x),len(y))-1))-1)
    x = x + [0] * (n - len(x))
    y = y + [0] * (n - len(y))
    z = convolve(x, y, b-1)
    # just carry forward to make sure no digit is more than b
    carry = 0
    non_zero = 0
    for i in range(len(z)):
        carry += z[i]
        z[i] = carry % b
        if z[i] != 0:
            non_zero = i
        carry = carry / b
    # carry should be zero here
    # shave off trailing zeros
    return z[:non_zero+1]

def multiply(x, y):
    X = [int(digit) for digit in bin(x)[:1:-1]]
    Y = [int(digit) for digit in bin(y)[:1:-1]]
    Z = multiply_digits(X, Y, 2)
    out = 0
    for bit in reversed(Z):
        out = (out << 1) | bit
    return out

assert multiply(13, 10) == 130

import random
random.seed(2478721839)

# generate random digit array of length n in base b
def generate(n, b = 2):
    return [random.randint(0, b-1) for i in range(n)]

# for digit arrays a, b, c: test associative property: a * (b * c) == (b * a) * c
def test(x, y, z, b = 2):
    return multiply_digits(x, multiply_digits(y, z, b), b) == multiply_digits(multiply_digits(y, x, b), z, b)

def test_size(n, b = 2):
    return test(generate(n, b), generate(n, b), generate(n, b), b)

import time

for i in range(20):
    t = time.time()
    result = test_size(1<<i)
    t = (time.time() - t)/4 # 4 multiplications done
    print 'Test 2^' + str(i) + ' digits | ' + \
          ('Success' if result else 'Failure') + ' | ' + "{0:.2f}".format(t) + ' seconds'

'''
Test 2^0 digits | Success | 0.00 seconds
Test 2^1 digits | Success | 0.00 seconds
Test 2^2 digits | Success | 0.00 seconds
Test 2^3 digits | Success | 0.00 seconds
Test 2^4 digits | Success | 0.00 seconds
Test 2^5 digits | Success | 0.00 seconds
Test 2^6 digits | Success | 0.01 seconds
Test 2^7 digits | Success | 0.01 seconds
Test 2^8 digits | Success | 0.02 seconds
Test 2^9 digits | Success | 0.04 seconds
Test 2^10 digits | Success | 0.06 seconds
Test 2^11 digits | Success | 0.17 seconds
Test 2^12 digits | Success | 0.34 seconds
Test 2^13 digits | Success | 0.76 seconds
Test 2^14 digits | Success | 1.40 seconds
Test 2^15 digits | Success | 3.95 seconds
Test 2^16 digits | Success | 8.52 seconds
Test 2^17 digits | Success | 18.07 seconds
Test 2^18 digits | Success | 37.90 seconds
Test 2^19 digits | Success | 80.93 seconds
'''

# Testing larger bases
for i in range(20):
    t = time.time()
    result = test_size(1<<10, 1<<(i*2))
    t = (time.time() - t)/4 # 4 multiplications done
    print 'Test base 2^' + str(i*2) + ' | ' + \
          ('Success' if result else 'Failure') + ' | ' + "{0:.2f}".format(t) + ' seconds'

'''
Test base 2^0 | Success | 0.05 seconds
Test base 2^2 | Success | 0.06 seconds
Test base 2^4 | Success | 0.09 seconds
Test base 2^6 | Success | 0.09 seconds
Test base 2^8 | Success | 0.12 seconds
Test base 2^10 | Success | 0.18 seconds
Test base 2^12 | Success | 0.17 seconds
Test base 2^14 | Success | 0.17 seconds
Test base 2^16 | Success | 0.20 seconds
Test base 2^18 | Success | 0.51 seconds
Test base 2^20 | Success | 0.26 seconds
Test base 2^22 | Success | 0.21 seconds
Test base 2^24 | Success | 0.29 seconds
Test base 2^26 | Success | 0.20 seconds
Test base 2^28 | Success | 0.23 seconds
Test base 2^30 | Success | 0.30 seconds
Test base 2^32 | Success | 0.45 seconds
Test base 2^34 | Success | 0.24 seconds
Test base 2^36 | Success | 0.72 seconds
Test base 2^38 | Success | 0.60 seconds
'''

t = time.time()
result = test_size(1<<15, 1<<(1<<5))
t = (time.time() - t)/4 # 4 multiplications done
print 'Test 2^15 groups of 2^5 digits' + ' | ' + \
      ('Success' if result else 'Failure') + ' | ' + "{0:.2f}".format(t) + ' seconds'
# 10.05 seconds instead of 90 seconds using no bit grouping
