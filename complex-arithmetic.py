# We use FFT on double precision complex numbers using complex n-th root of unity

from cmath import *

def unity(n):
    return complex(cos(2*pi/n), sin(2*pi/n))

# forward fft using cooley tukey algorithm
def fft(sequence, root):
    if (len(sequence) == 1):
        return
    even = sequence[::2]
    odd = sequence[1::2]
    r = root * root
    fft(even, r)
    fft(odd, r)
    r = 1
    m = len(sequence)
    for i in range(m):
        j = i % (m/2)
        sequence[i] = even[j] + odd[j] * r
        r *= root

# inverse fft
def ifft(sequence, root):
    m = len(sequence)
    for i in range(1, (m+1)/2):
        sequence[m-i], sequence[i] = sequence[i], sequence[m-i]
    fft(sequence, root)
    for i in range(m):
        sequence[i] /= m
        sequence[i] = int(sequence[i].real+0.5)

def convolve(x, y):
    n = len(x)
    r = unity(n)
    fft(x, r)
    fft(y, r)
    z = [x[i] * y[i] for i in range(len(x))]
    ifft(z, r)
    return z

# x, y are digit arrays (low to high order) in base b
# X = x[0] + x[1]b + x[1]b^2 + ...
def multiply_digits(x, y, b = 2):
    # round up digits to power of 2 with the later half all zeros
    n = 1 << (len(bin(max(len(x),len(y))-1))-1)
    x = x + [0] * (n - len(x))
    y = y + [0] * (n - len(y))
    z = convolve(x, y)
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
Test 2^6 digits | Success | 0.00 seconds
Test 2^7 digits | Success | 0.01 seconds
Test 2^8 digits | Success | 0.01 seconds
Test 2^9 digits | Success | 0.03 seconds
Test 2^10 digits | Success | 0.05 seconds
Test 2^11 digits | Success | 0.10 seconds
Test 2^12 digits | Success | 0.23 seconds
Test 2^13 digits | Success | 0.50 seconds
Test 2^14 digits | Success | 1.07 seconds
Test 2^15 digits | Success | 2.16 seconds
Test 2^16 digits | Success | 4.66 seconds
Test 2^17 digits | Success | 9.76 seconds
Test 2^18 digits | Success | 20.80 seconds
Test 2^19 digits | Success | 44.68 seconds
'''

# Fails due to numerical precision reasons for larger bases
print test_size(1<<10, 1<<15) # True
print test_size(1<<10, 1<<16) # False
print test_size(1<<15, 1<<11) # True
print test_size(1<<15, 1<<12) # False
