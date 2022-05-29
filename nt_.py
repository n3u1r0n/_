try:
  import cupy as lib
except:
  import numpy as lib
import numpy as np
import os


Primes = lib.array([], dtype=int)
SmallestFactors = lib.array([], dtype=int)


def sieve_primes(n):
  global Primes
  if n <= 3: Primes = lib.array([2, 3], dtype=int)[:max(0, n - 1)]
  else:
    n += 1
    sieve_array = lib.ones(n // 3 + (n % 6 == 2), dtype=bool)
    sieve_array[0] = 0
    i = 0
    while 9 * i * i <= n:
      if sieve_array[i]:
        k = (3 * i + 1) | 1
        sieve_array[((k * k) // 3)::2 * k] = 0
        sieve_array[(k * k + 4 * k - 2 * k * (i & 1)) // 3::2 * k] = 0
      i += 1
    Primes = lib.r_[2, 3, 3 * lib.flatnonzero(sieve_array) + 1 | 1]
    del sieve_array

def sieve_smallest_factors(n):
  global SmallestFactors
  if n <= 3: SmallestFactors = lib.array([0, 1, 2, 3], dtype=int)[:max(0, n + 1)]
  else:
    sieve_array = lib.arange(n + 1, dtype=int)
    sieve_array[2::2] = 2
    sieve_array[3::6] = 3
    i = 1
    while (k := i * 3 + 1 | 1) * k <= n:
      sieve_array[k::k] = lib.minimum(sieve_array[k::k], k)
      i += 1
    SmallestFactors = sieve_array
    del sieve_array

def sieve_all(n):
  sieve_primes(n)
  sieve_smallest_factors(n)

def load_primes(filename='primes.npy'):
  global Primes
  assert type(filename) == str and os.path.isfile(filename)
  Primes = lib.load(filename, allow_pickle=True)

def save_primes(filename='primes.npy'):
  assert type(filename) == str
  lib.save(filename, Primes)

def load_smallest_factors(filename='smallest_factors.npy'):
  global SmallestFactors
  assert type(filename) == str and os.path.isfile(filename)
  SmallestFactors = lib.load(filename, allow_pickle=True)

def save_smallest_factors(filename='smallest_factors.npy'):
  assert type(filename) == str
  lib.save(filename, SmallestFactors)

def pi(x):
  assert 0 <= x <= Primes[-1]
  if x < 2: return 0
  low, high = 0, len(Primes)
  while high - low > 1:
    center = (high + low) // 2
    if x < Primes[center]: high = center
    else: low = center
  return low + 1

def is_prime(n):
  assert 0 <= n <= Primes[-1]
  if n == 2 or n == 3: return True
  if n < 2 or n % 2 == 0 or n % 3 == 0: return False
  return bool(Primes[pi(n) - 1] == n)

def gcd(a, b):
  while b != 0:
    a, b = b, a - b * (a // b)
  return abs(a)

def prime_factor_decomposition(n):
  assert n < len(SmallestFactors)
  prime_factors = []
  while n != 1:
    prime_factors += [k := SmallestFactors[n]]
    n //= k
  return lib.array(prime_factors, dtype=int)

def mod_inverse(n, m):
  return pow(n, -1, m)

def sqrt(n):
  assert type(n) == int and n >= 0
  a = n
  while a * a > n:
    a = (a + n // a) // 2
  return a

def is_square(n):
  if type(n) != int or n < 0: return False
  return (root := sqrt(n)) * root == n