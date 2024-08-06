import math
import numpy  
import sympy as sym
from sympy import Lambda
from sympy import Pow
from sympy import arg
import time


print("coucou")


f = sym.Symbol('f')
x = sym.Symbol('x')

deg = 1
deg = deg * (1-(x* -0.4)) + f*f
deg = deg * (1-(x* -0.4)) + f
deg = deg * (1-(x* -0.4))

deg = sym.expand(deg)

print(deg)

#deg = deg.replace(lambda a: Pow(x,a), Pow(x,1))

#deg = deg.replace(lambda x: x.is_Mul, lambda x: 2*x)
#deg = deg.replace(x.is_Pow, Pow(x,1))


#deg = deg.replace(lambda a: Pow(x,a), lambda a: Pow(x,1))
deg = deg.replace(lambda a: a.is_Pow and a.base == x, lambda a: Pow(x,1))

#deg = deg.replace(Pow, lambda a,b: Pow(a,1))

print(deg)