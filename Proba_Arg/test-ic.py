import math
import numpy  
import sympy as sym
from sympy import Lambda
from sympy import Pow
from sympy import arg
import time

import pickle

#from sage.all import *
#from sage.arith.power import generic_power

import matplotlib.pyplot as plt

"""
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
"""

som = pickle.load(open("TEST2_max_50_70_1500", "rb"))
#temps = pickle.load(open("TEST2_max_50_70_1500", "rb"))

print(som)
#print(temps)

print(len(som))
#print(len(temps))

for i in range(1,100):
    som.append(i)

pickle.dump(som,open("TEST2_max_50_70_1500", "wb"))
#pickle.dump(open("TEST2_max_50_70_1500", "wb"))

som = pickle.load(open("TEST2_max_50_70_1500", "rb"))
print(som)
print(len(som))

#plt.plot(som,temps,"ob") # ob = type de points "o" ronds, "b" bleus
#plt.ylabel('Times')
#plt.xlabel('Sum Termes')
#plt.show()

