from helper import trim_zeros, fill_zeros
from functools import reduce
from nt import gcd

class Polynomial:
  def __init__(self, *coeffs):
    self.coeffs = trim_zeros(list(coeffs))
    if len(self.coeffs) == 0: self.coeffs = [0]
    self.degree = len(self.coeffs) - 1

  def __str__(self):
    if self.degree == 0: return str(self.coeffs[0])
    string = ''
    for n, coeff in enumerate(self.coeffs):
      if coeff == 0: pass
      elif n == 0:
        if coeff == 1: string += '1'
        elif coeff == -1: string += '-1'
        else: string += str(coeff)
      elif (not isinstance(coeff, complex)) and coeff < 0:
        string += '-X' if coeff == -1 else (str(coeff) + 'X')
        string += ('^{' + str(n) + '}') if n > 1 else ''
      else:
        string += '+' if string != '' else ''
        string += (str(coeff) + 'X') if coeff != 1 else 'X'
        string += ('^{' + str(n) + '}') if n > 1 else ''
    return string.replace('j', 'i')

  def __add__(self, other):
    if isinstance(other, Polynomial):
      coeffs = [
        self_coeff + other_coeff for self_coeff, other_coeff in
        zip(fill_zeros(self.coeffs, other.degree + 1), fill_zeros(other.coeffs, self.degree + 1))
      ]
    else:
      coeffs = [self.coeffs[0] + other] + self.coeffs[1:]
    return Polynomial(*coeffs)
  
  __radd__ = __add__

  def __sub__(self, other):
    if isinstance(other, Polynomial):
      coeffs = [
        self_coeff - other_coeff for self_coeff, other_coeff in
        zip(fill_zeros(self.coeffs, other.degree + 1), fill_zeros(other.coeffs, self.degree + 1))
      ]
    else:
      coeffs = [self.coeffs[0] - other] + self.coeffs[1:]
    return Polynomial(*coeffs)
  
  def __neg__(self):
    coeffs = [- coeff for coeff in self.coeffs]
    return Polynomial(*coeffs)
  
  def __rsub__(self, other):
    if isinstance(other, Polynomial):
      coeffs = [
        other_coeff - self_coeff for self_coeff, other_coeff in
        zip(fill_zeros(self.coeffs, other.degree + 1), fill_zeros(other.coeffs, self.degree + 1))
      ]
    else:
      coeffs = [other - self.coeffs[0]] + [-coeff for coeff in self.coeffs[1:]]
    return Polynomial(*coeffs)

  def __mul__(self, other):
    if isinstance(other, Polynomial):
      coeffs = [
        sum([
          self.coeffs[m] * other.coeffs[n - m] 
          for m in range(max(0, n - other.degree), min(n, self.degree) + 1)
        ]) for n in range(self.degree + other.degree + 1)
      ]
    else:
      coeffs = [other * coeff for coeff in self.coeffs]
    return Polynomial(*coeffs)
  
  __rmul__ = __mul__

  __repr__ = __str__

  def __eq__(self, other):
    if isinstance(other, Polynomial):
      return self.degree == other.degree and all([
        self_coeff == other_coeff for self_coeff, other_coeff in zip(self.coeffs, other.coeffs)
      ])
    return self.degree == 0 and self.coeffs[0] == other
  
  def __divmod__(self, other):
    if isinstance(other, Polynomial):
      remainder = self.coeffs[:]
      coeffs = []
      for n in reversed(range(self.degree - other.degree + 1)):
        coeffs = [remainder[-1] / other.coeffs[-1]] + coeffs
        remainder = remainder[:n] + [remainder[i + n] - coeffs[0] * other.coeffs[i] for i in range(other.degree)]
      return Polynomial(*coeffs), Polynomial(*remainder)
    coeffs = [coeff / other for coeff in self.coeffs]
    return Polynomial(*coeffs), Polynomial(0)
  
  def __floordiv__(self, other):
    if isinstance(other, Polynomial):
      remainder = self.coeffs[:]
      coeffs = []
      for n in reversed(range(self.degree - other.degree + 1)):
        coeffs = [remainder[-1] / other.coeffs[-1]] + coeffs
        remainder = remainder[:n] + [remainder[i + n] - coeffs[0] * other.coeffs[i] for i in range(other.degree)]
      return Polynomial(*coeffs)
    coeffs = [coeff / other for coeff in self.coeffs]
    return Polynomial(*coeffs)
  
  def __mod__(self, other):
    if isinstance(other, Polynomial):
      remainder = self.coeffs[:]
      coeffs = []
      for n in reversed(range(self.degree - other.degree + 1)):
        coeffs = [remainder[-1] / other.coeffs[-1]] + coeffs
        remainder = remainder[:n] + [remainder[i + n] - coeffs[0] * other.coeffs[i] for i in range(other.degree)]
      return Polynomial(*remainder)
    coeffs = [coeff / other for coeff in self.coeffs]
    return Polynomial(0)

  def __truediv__(self, other):
    if isinstance(other, Polynomial):
      return Rational(self, other)
    coeffs = [coeff / other for coeff in self.coeffs]
    return Polynomial(*coeffs)

  def __call__(self, x):
    return reduce(lambda y, z: x * y + z, self.coeffs[::-1])
  
  def __abs__(self):
    return Polynomial(*self.coeffs)
  
  def __pow__(self, n):
    if n == 0:
      return Polynomial(1)
    return self * self ** (n - 1)

  @property
  def prime(self):
    return Polynomial(*[(n + 1) * coeff for n, coeff in enumerate(self.coeffs[1:])])
  
  @property
  def leading(self):
    return self.coeffs[-1]
  
  @property
  def constant(self):
    return self.coeffs[0]
  
  @staticmethod
  def Monomial(n):
    return Polynomial(*([0] * n), 1)

X = Polynomial(0, 1)

class Rational:
  def __init__(self, num, den):
    if (GCD := gcd(num, den)).degree > 0:
      self.num = num // GCD
      self.den = den // GCD
    else:
      self.num = num
      self.den = den
  
  def __str__(self):
    return '\\frac{' + str(self.num) + '}{' + str(self.den) + '}'
