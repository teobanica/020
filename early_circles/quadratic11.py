from cyclotomic60 import K,Q,POW
class L:
 def __init__(self,a=0,b=0):
  if isinstance(a,L):self.a,self.b=a.a,a.b;return
  self.a,self.b=K(a),K(b)
 def __add__(x,y):y=L(y);return L(x.a+y.a,x.b+y.b)
 __radd__=__add__
 def __neg__(x):return L(-x.a,-x.b)
 def __sub__(x,y):return x+-L(y)
 def __rsub__(x,y):return L(y)+-x
 def __mul__(x,y):y=L(y);return L(x.a*y.a+11*x.b*y.b,x.a*y.b+x.b*y.a)
 __rmul__=__mul__
 def __bool__(x):return bool(x.a) or bool(x.b)
 def __eq__(x,y):y=L(y);return x.a==y.a and x.b==y.b
 def inverse(x):
  n=x.a*x.a-11*x.b*x.b;assert n
  return L(x.a/n,-x.b/n)
 def __truediv__(x,y):return x*L(y).inverse()
 def conj(x):return L(x.a.conj(),x.b.conj())
 def serial(x):return [x.a.serial(),x.b.serial()]
def rrefL(A):
 a=[[L(x) for x in row] for row in A];piv=[];r=0
 for c in range(len(a[0])):
  s=next((i for i in range(r,len(a)) if a[i][c]),None)
  if s is None:continue
  a[r],a[s]=a[s],a[r];v=a[r][c].inverse();a[r]=[x*v for x in a[r]]
  for i in range(len(a)):
   if i!=r and a[i][c]:
    v=a[i][c];a[i]=[x-v*y for x,y in zip(a[i],a[r])]
  piv.append(c);r+=1
  if r==len(a):break
 return a,piv
