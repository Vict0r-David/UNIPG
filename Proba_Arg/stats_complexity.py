import math
import numpy  
import sympy as sym
from sympy import Lambda
from sympy import Pow
from sympy import arg
import time

import pickle

from sage.all import *
from sage.arith.power import generic_power

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

#max1 = pickle.load(open("save_max_50_70_1500", "rb"))
#time1 = pickle.load(open("save_time_50_70_1500", "rb"))

#max2 = pickle.load(open("max_50_75_1000", "rb"))
#time2 = pickle.load(open("time_50_75_1000", "rb"))

#max12 = pickle.load(open("merge_max_50_70-5_15-1000","rb"))
#time12 = pickle.load(open("merge_time_50_70-5_15-1000","rb"))

maxi = pickle.load(open("max123_50_70-5_15-1000", "rb"))
time = pickle.load(open("time123_50_70-5_15-1000", "rb"))

ind = maxi.index(559872)

del maxi[ind]
del time[ind]

#maxi = max1 + max2
#time = time1 + time2

#pickle.dump(maxi,open("merge_max_50_70-5_15-1000", "wb"))
#pickle.dump(time,open("merge_time_50_70-5_15-1000", "wb"))

#maxi = pickle.load(open("merge_max_50_70-5_15-1000", "rb"))
#time = pickle.load(open("merge_time_50_70-5_15-1000", "rb"))

#maxi = max12 + max3
#time = time12 + time3

pickle.dump(maxi,open("finalMax_50_70-5_15-1000", "wb"))
pickle.dump(time,open("finalTime_50_70-5_15-1000", "wb"))

print(f"nombre de valeurs {len(time)}")

dico = {}
dico[1] = 0
dico[10] = 0
dico[20] = 0
dico[30] = 0
dico[40] = 0
dico[50] = 0
dico[60] = 0
dico[70] = 0
dico[80] = 0
dico[90] = 0

for i in range(0,10700,100):
    dico[i] = 0

print(dico)

for elem in time:
    #print(int(elem))
    if int(elem) < 1:
        elem = 0
    elif int(elem) < 10:
        elem = 1
    elif int(elem) < 20:
        elem = 10
    elif int(elem) < 30:
        elem = 20
    elif int(elem) < 40:
        elem = 30
    elif int(elem) < 50:
        elem = 40
    elif int(elem) < 60:
        elem = 50
    elif int(elem) < 70:
        elem = 60
    elif int(elem) < 80:
        elem = 70
    elif int(elem) < 90:
        elem = 80
    elif int(elem) < 100:
        elem = 90
    else:
        modu = int(elem) % 100
        elem = int(elem) - modu

    dico[int(elem)] += 1

liste_nb_time = []

for key,elem in dico.items():
    liste_nb_time.append([elem/len(time),key]) 

print(dico)

print(liste_nb_time)

"""
{0: 124268, 100: 3, 200: 2, 300: 0, 400: 2, 500: 0, 600: 0, 700: 0, 800: 1, 900: 1, 1000: 0, 1100: 1, 1200: 0, 1300: 0, 1400: 1, 1500: 0, 1600: 0, 1700: 0, 1800: 0, 1900: 0, 2000: 0, 2100: 0, 2200: 0, 2300: 0, 2400: 0, 
2500: 1, 2600: 0, 2700: 0, 2800: 0, 2900: 0, 3000: 0, 3100: 0, 3200: 0, 3300: 0, 3400: 0, 3500: 0, 3600: 0, 3700: 0, 3800: 0, 3900: 0, 4000: 0, 4100: 0, 4200: 0, 4300: 0, 4400: 0, 4500: 0, 4600: 0, 4700: 0, 4800: 0, 4900: 0, 5000: 1, 5100: 0}
[[0.9998953983312011, 0], [2.413884664590726e-05, 100], [1.609256443060484e-05, 200], [0.0, 300], [1.609256443060484e-05, 400], [0.0, 500], [0.0, 600], [0.0, 700], [8.04628221530242e-06, 800], [8.04628221530242e-06, 900], [0.0, 1000], [8.04628221530242e-06, 1100], [0.0, 1200], [0.0, 1300], [8.04628221530242e-06, 1400], [0.0, 1500], [0.0, 1600], [0.0, 1700], [0.0, 1800], [0.0, 1900], [0.0, 2000], [0.0, 2100], [0.0, 2200], [0.0, 2300], [0.0, 2400], [8.04628221530242e-06, 2500], [0.0, 2600], [0.0, 2700], [0.0, 2800], [0.0, 2900], [0.0, 3000], [0.0, 3100], [0.0, 3200], [0.0, 3300], [0.0, 3400], [0.0, 3500], [0.0, 3600], [0.0, 3700], [0.0, 3800], [0.0, 3900], [0.0, 4000], [0.0, 4100], [0.0, 4200], [0.0, 4300], [0.0, 4400], [0.0, 4500], [0.0, 4600], [0.0, 4700], [0.0, 4800], [0.0, 4900], [8.04628221530242e-06, 5000], [0.0, 5100]]

nombre de valeurs 124281
 0: 124235, 10: 11, 20: 9, 30: 3, 40: 1, 50: 2, 60: 2, 70: 1, 80: 2, 90: 2, 100: 3, 
 200: 2, 300: 0, 400: 2, 500: 0, 600: 0, 700: 0, 800: 1, 900: 1, 1000: 0, 1100: 1,
 1400 : 1, 2500 : 1, 5000: 1
"""

plt.plot(maxi,time,"ob") # ob = type de points "o" ronds, "b" bleus
plt.ylabel('Times (seconds)')
plt.xlabel('Maximum number of terms among all dependent argument resolutions')
plt.show()

#plt.plot(som,temps,"ob") # ob = type de points "o" ronds, "b" bleus
#plt.ylabel('Times')
#plt.xlabel('Sum Termes')
#plt.show()


"""

nombre de valeurs 198872
 0: 198791, 10: 16, 20: 12, 30: 4, 40: 3, 50: 4, 60: 3, 70: 3, 80: 2, 90: 6, 100: 10, 
 200: 3, 300: 0, 400: 2, 500: 1, 600: 0, 700: 0, 800: 2, 900: 1, 1000: 1, 1100: 1,
 1400: 1, 2500: 1, 4400: 1, 5000: 1 , 10600: 1
"""

