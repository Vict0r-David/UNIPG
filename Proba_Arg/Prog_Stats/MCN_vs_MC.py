
# Created 03/10/2022

#import math
import numpy  

import time

from sage.all import *
from sage.arith.power import generic_power

import matplotlib.pyplot as plt

import pickle

import MonteCarlo as MC


List_diff = pickle.load(open("diff_50_75_1000_rapide", "rb"))
List_error = pickle.load(open("error_50_75_1000_rapide", "rb"))
l_error_small = pickle.load(open("error_borne_50_75_1000_rapide", "rb"))



figure = plt.figure(figsize = (10, 10))
plt.gcf().subplots_adjust(left = 0.2, bottom = 0.2,
                       right = 0.7, top = 0.7, 
                       wspace = 0.5, hspace = 0.7)

axes = figure.add_subplot(1, 2, 1)
plt.hist(List_diff,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
#plt.title("Histogram of 1000 MCN (50 arg, 75 att)")
plt.xlabel('Error difference')

#plt.yticks([0,1000,3000,5000,10000,15000,20000,25000,30000])
plt.ylabel('Number of occurrences')


axes = figure.add_subplot(1, 2, 2)
#plt.hist(List_error,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
plt.hist(l_error_small,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
#plt.title("Histogram of 1000 MCN (50 arg, 75 att)")

#plt.yticks([0,1000,3000,5000,10000,15000,20000])
plt.xlabel('Error percentage')
#plt.ylabel('Number of occurrence')

plt.show()


"""
avg diff = 0.16035736307631446
avg error = 39.391539887362455
approximation == 0, = 12139
proba == approximation, = 13450
cptProb_0 et approx not 0, = 0
nb valeur error percentage = 37582
"""


############################### COMPLEXITY
"""
plt.plot(tab_Max,tab_Time,"ob") # ob = type de points "o" ronds, "b" bleus
plt.ylabel('Times')
plt.xlabel('Max Termes')
plt.show()


plt.plot(tab_Sum,tab_Time,"ob") # ob = type de points "o" ronds, "b" bleus
plt.ylabel('Times')
plt.xlabel('Sum Termes')
plt.show()
"""

