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
file_AF = ".\SCN_500\SCN_500_1.txt"
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
"""
# Quicksort Sort

# This implementation utilizes pivot as the last element in the nums list
# It has a pointer to keep track of the elements smaller than the pivot
# At the very end of partition() function, the pointer is swapped with the pivot
# to come up with a "sorted" nums relative to the pivot

# Function to find the partition position
def partition(array, low, high):
    # choose the rightmost element as pivot
    pivot = array[high][1]
    # pointer for greater element
    i = low - 1
    # traverse through all elements
    # compare each element with pivot
    for j in range(low, high):
        if array[j][1] <= pivot:
            # If element smaller than pivot is found
            # swap it with the greater element pointed by i
            i = i + 1
            # Swapping element at i with element at j
            (array[i], array[j]) = (array[j], array[i])
    # Swap the pivot element with the greater element specified by i
    (array[i + 1], array[high]) = (array[high], array[i + 1])
    # Return the position from where partition is done
    return i + 1

# function to perform quicksort
def quickSort(array, low, high):
    if low < high:
        # Find pivot element such that
        # element smaller than pivot are on the left
        # element greater than pivot are on the right
        pi = partition(array, low, high)
        # Recursive call on the left of pivot
        quickSort(array, low, pi - 1)
        # Recursive call on the right of pivot
        quickSort(array, pi + 1, high)
"""

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
        level_Arg(att, lvl+1, ldico_att_in,dico_lvl)
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


"""
def update(dico_lvl,b,val,s_att):
    diff = val - dico_lvl[b]
    dico_lvl[b] = val
    for a in dico_lvl:
        if b in dico_lvl[a][1]:
            dico_lvl[b][0] += diff
    return dico_lvl

def level(arg,d_att):
    s_att = set(d_att[arg])
    dico_lvl = {}
    dico_lvl[arg] = [0, []]
    i = 0
    while i < len(s_att):
        for b in d_att[s_att[i]]:
            if b not in dico_lvl:
                dico_lvl[b] = [dico_lvl[s_att[i]][0] + 1, dico_lvl[s_att[i]][1].add(s_att[i]) ]
            elif dico_lvl[b] < dico_lvl[s_att[i]][0] + 1:
                dico_lvl = update(dico_lvl, b, dico_lvl[s_att[i]][0] + 1, s_att)
        i += 1
    output = {}
    for a,info in dico_lvl:
        output[a] = info[0]
    heapSort(output)
    return output


def max_len(a,list_Arg,list_Att):
    G = nx.DiGraph()
    G.add_nodes_from(list_Arg)
    G.add_edges_from(list_Att)
    list_lvl = [[a,0]] 
    #print(list_Arg)
    list_Arg.remove(a)
    #print(list_Arg)
    for b in list_Arg:
        #print(b)
        #print(list_lvl)
        try:
            length = len(max(nx.all_simple_paths(G, a, b), key=lambda x: len(x)))
            list_lvl.append(b,length)
        except:
            continue
    return list_lvl

#print(max_len("a203",list_Arg,list_Att))
"""

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

#dico_lvl = init_dico_lvl(dico_Arg,dico_lvl)
start = time.time()
proba = fastSCN("a1",dico_Arg,dico_Att, dico_lvl)
end = time.time()
elapsed = end - start

print(f"Probability of a (SCN algo) = {proba}")
print(f"Time (SCN algo) = {elapsed}\n")