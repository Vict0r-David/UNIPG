import math
import numpy  
import sympy as sym
from sympy import Lambda
from sympy import Pow
from sympy import arg
import time

import networkx as nx
import random
from itertools import groupby


import sys
 
# the setrecursionlimit function is
# used to modify the default recursion
# limit set by python. Using this,
# we can increase the recursion limit
# to satisfy our needs
 
sys.setrecursionlimit(10**6)

#### GESTION DU GRAPHE  - INITIALISATION ####

#file_AF = "AF_test2.txt"
#file_AF = ".\DAG_45\DAG_45_0.06_1.txt"
file_AF = ".\SCN_800\SCN_800_1.txt"
#file_AF = ".\Old_test\AF5_3.txt"
AF = open(file_AF,"r")
dico_Arg = {}
dico_Att = {}
dico_lvl = {}
list_Arg = []
list_Att = []

for line in AF:
    if line[0:3] == "arg":
        l1 = line.partition("(")
        l2 = l1[2].partition(")")
        list_Arg.append(l2[0])
        dico_Arg[l2[0]] = 1
        dico_lvl[l2[0]] = 0

    if line[0:3] == "att":
        l1 = line.partition("(")
        l2 = l1[2].partition(")")
        l3 = l2[0].partition(",")
        att_from = l3[0]
        att_to = l3[2]
        l4 = l2[2][1:].partition("\n")
        #w = float(l4[0][:-1])
        w = numpy.float64(l4[0][:-1])
        att = [att_from, att_to, w]
        list_Att.append((att_from,att_to))
        id = att_from+"->"+att_to
        dico_Att[id] = att

#print(dico_Arg)
#print(dico_Att)

#AF_1.txt
#dico_Arg = {"a":1, "b":1, "c":1, "d":1}
#dico_Att = {"a->b":["a","b",-0.6], "b->c": ["b","c",-0.3], "d->b":["d","b",-0.2], "a->d":["a","d",-0.5], "c->a":["c","a",-0.2]}

#AF_2.txt
#dico_Arg = {"a":1, "b":1, "c":1, "d":1, "e":1}
#dico_Att = {"a->c": ["a","c",-0.3], "b->c": ["b","c",-0.9], "c->e": ["c","e",-0.4], "d->e": ["d","e",-0.3]}



def list_dico_Att_in(dico_Arg, dico_Att):
    ldico_att_in = {}
    for id_arg in dico_Arg:
        ldico_att_in[id_arg] = []
    for id,att in dico_Att.items():
        ldico_att_in[att[1]] += [id]
    return ldico_att_in


def dlist_Att_to_Arg(dico_list_att):
    dico = {}
    for arg,l in dico_list_att.items():
        dico[arg] = []
        for att in l:
            l2 = att.split("-")
            dico[arg].append(l2[0])
    return dico


################################################################################################################
################################################################################################################

############# Heap sort


# To heapify subtree rooted at index i.
# n is size of heap
  
  
def heapify(arr, n, i):
    largest = i  # Initialize largest as root
    l = 2 * i + 1  # left = 2*i + 1
    r = 2 * i + 2  # right = 2*i + 2
    # See if left child of root exists and is
    # greater than root
    if l < n and arr[i][1] < arr[l][1]:
        largest = l
  
    # See if right child of root exists and is
    # greater than root
    if r < n and arr[largest][1] < arr[r][1]:
        largest = r
  
    # Change root, if needed
    if largest != i:
        (arr[i], arr[largest]) = (arr[largest], arr[i])  # swap
        # Heapify the root.
        heapify(arr, n, largest)
  
  
# The main function to sort an array of given size
def heapSort(arr):
    n = len(arr)
    # Build a maxheap.
    # Since last parent will be at ((n//2)-1) we can start at that location.
    for i in range(n // 2 - 1, -1, -1):
        heapify(arr, n, i)
    # One by one extract elements
    for i in range(n - 1, 0, -1):
        (arr[i], arr[0]) = (arr[0], arr[i])  # swap
        heapify(arr, i, 0)
  
#arr = [["a",5],["a",10],["a",2],["a",0]]
#heapSort(arr)
#print(arr)

############################################### Acyclic Case ##############################################

# Preprocessing

def level_Arg(start, lvl, ldico_att_in,dico_lvl):
    for att in ldico_att_in[start]:
        if dico_lvl[att] < lvl:
            dico_lvl[att] = lvl
        level_Arg(att, lvl+1, ldico_att_in, dico_lvl)
    return dico_lvl


def init_dico_lvl(dico_Arg,dico_lvl):
    for arg in dico_Arg:
        dico_lvl[arg] = 0
    return dico_lvl


def order_level(start,dico_lvl):
    liste = []
    for arg,val in dico_lvl.items():
            if val > 0 or arg == start:
                liste.append([arg,val])
    #size = len(liste)
    #quickSort(liste, 0, size - 1)
    heapSort(liste)
    return liste


def fastSCN(goal,dico_Arg,dico_Att,dico_lvl):
    ldico_att_in = list_dico_Att_in(dico_Arg, dico_Att)
    dico_att_in = dlist_Att_to_Arg(ldico_att_in)
    dico_lvl = level_Arg(goal, 1, dico_att_in, dico_lvl)
    liste_lvl = order_level(goal,dico_lvl)
    #l_lvl = level(goal,dico_att_in)

    dico_deg = {}
    i = len(liste_lvl)-1
    
    while i >= 0: 
        l_att = ldico_att_in[str(liste_lvl[i][0])]
        if l_att == []:
            dico_deg[str(liste_lvl[i][0])] = 1
        else:
            deg = 1
            for att in l_att:
                deg = deg * (1-(dico_deg[dico_Att[att][0]]*- dico_Att[att][2]))
            dico_deg[str(liste_lvl[i][0])] = deg
        i-=1
    return dico_deg[goal]

"""
#dico_lvl = init_dico_lvl(dico_Arg,dico_lvl)
start = time.time()
proba = fastSCN("a1",dico_Arg,dico_Att, dico_lvl)
end = time.time()
elapsed = end - start

print(f"Probability of a (SCN algo) = {proba}")
print(f"Time (SCN algo) = {elapsed}\n")
"""