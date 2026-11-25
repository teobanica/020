"""Exact rational interval arithmetic; no floating-point operations."""
from fractions import Fraction as Q
from math import isqrt
class I:
 def __init__(self,a,b=None):self.a=Q(a);self.b=Q(a if b is None else b);assert self.a<=self.b
 def __add__(x,y):
  if not isinstance(y,I):y=I(y)
  return I(x.a+y.a,x.b+y.b)
 __radd__=__add__
 def __neg__(x):return I(-x.b,-x.a)
 def __sub__(x,y):return x+-asI(y)
 def __rsub__(x,y):return asI(y)+-x
 def __mul__(x,y):
  y=asI(y);v=[x.a*y.a,x.a*y.b,x.b*y.a,x.b*y.b];return I(min(v),max(v))
 __rmul__=__mul__
 def sqrt(x,digits):
  assert x.a>=0;n=10**digits
  l=isqrt(x.a.numerator*n*n//x.a.denominator)
  h=isqrt(x.b.numerator*n*n//x.b.denominator)+1
  assert Q(l*l,n*n)<=x.a and Q(h*h,n*n)>x.b
  return I(Q(l,n),Q(h,n))
 def serial(x):return [str(x.a),str(x.b)]
def asI(x):return x if isinstance(x,I) else I(x)
def real_bound(coeffs,digits=4):
 s5=I(5).sqrt(digits);s3=I(3).sqrt(digits);s10=(10-2*s5).sqrt(digits)
 # zeta_60 = exp(i*(pi/5-pi/6)).
 re=((1+s5)*s3+s10)*Q(1,8)
 im=(s10*s3-(1+s5))*Q(1,8)
 a=I(0);b=I(0)
 for c in reversed(coeffs):a,b=a*re-b*im+Q(c),a*im+b*re
 return a
