# Victor DAVID
# Created 03/10/2022

import math
import numpy  
import sympy as sym
from sympy import Pow
import time


#### GESTION DU GRAPHE  - INITIALISATION ####

file_AF = "AF_0.txt"
AF = open(file_AF,"r")
dico_Arg = {}
dico_Att = {}
dico_lvl = {}
dico_pw = {}

for line in AF:
    if line[0:3] == "arg":
        l1 = line.partition("(")
        l2 = l1[2].partition(")")
        dico_Arg[l2[0]] = 1
        dico_lvl[l2[0]] = 0
        dico_pw[l2[0]] = 0
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



################################################################################################################
############################################### Usefull Function ###############################################



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



#############################################################################################
################################ INUTILE dans cette version #################################

def list_dico_Att_out(dico_Arg, dico_Att):
    ldico_att_out = {}
    for id_arg in dico_Arg:
        ldico_att_out[id_arg] = []
    for id,att in dico_Att.items():
        ldico_att_out[att[0]] += [id]
    return ldico_att_out


def get_paths(node, end, ldico_att_out, dico_att, paths=None, current_path=None):
    if paths is None:
        paths = []
    if current_path is None:
        current_path = []
    current_path.append(node)
    if node not in ldico_att_out or node == end:
        paths.append(current_path)
    else:
        children = []
        for att in ldico_att_out[node]:
            children.append(dico_att[att][1])
        for child in children:
            get_paths(child, end, ldico_att_out, dico_att, paths, list(current_path))
    return paths


def get_paths_back(node, ldico_att_in, dico_att, paths=None, current_path=None):
    if paths is None:
        paths = []
    if current_path is None:
        current_path = []
    current_path.append(node)
    if ldico_att_in[node] == []: 
        paths.append(current_path)
    else:
        children = []
        for att in ldico_att_in[node]:
            children.append(dico_att[att][0])
        for child in children:
            get_paths_back(child, ldico_att_in, dico_att, paths, list(current_path))
    return paths


def arg_to_att(path):
    output = []
    for i in range(len(path)):
        if i+1 < len(path):
            att = path[i] + '->' + path[i+1]
            output.append(att)
    return output


def good_paths(node, end, ldico_att_out, dico_att):
    all = get_paths(node, end, ldico_att_out, dico_att)
    paths = []
    for path in all:
        if path[-1] == end:
            p = arg_to_att(path)
            paths.append(p)
    return paths    


def all_att(node, end, ldico_att_out, dico_att):
    paths = good_paths(node, end, ldico_att_out, dico_att)
    set_att = set()
    for path in paths:
        for att in path:
            set_att.add(att)
    return set_att


def path_att_def(node, end, ldico_att_out, dico_att):
    paths = good_paths(node, end, ldico_att_out, dico_att)
    l_att = []
    l_def = []
    for path in paths:
        if len(path)%2 == 1:
            l_att.append(path)
        else:
            l_def.append(path)
    return l_att,l_def




################################################################################################################
################################################################################################################

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


############################################### Acyclic Case ##############################################

# Preprocessing

def level_Arg(start, lvl, ldico_att_in,dico_lvl):
    for att in ldico_att_in[start]:
        if dico_lvl[att] < lvl:
            dico_lvl[att] = lvl
        level_Arg(att, lvl+1, ldico_att_in,dico_lvl)
    return dico_lvl


def order_level(dico_lvl):
    liste = []
    for arg,val in dico_lvl.items():
        liste.append([arg,val])
    size = len(liste)
    quickSort(liste, 0, size - 1)
    return liste


# Symbolic Computation 

def develop(formula, arg, dico_att_in, dico_Att, dico_pw):
    new_arg = 1
    for att in dico_att_in[arg]:
        new_arg = new_arg * (1- sym.Symbol(dico_Att[att][0]) * float(-dico_Att[att][2]))
        dico_pw[dico_Att[att][0]] += 1
    dico_pw[arg] = 0
    formula = formula.subs(arg, new_arg)
    expF = sym.expand(formula)
    for att in dico_att_in[arg]:
        if dico_pw[dico_Att[att][0]] > 1:   
            for i in range(2,dico_pw[dico_Att[att][0]]+1):
                s = sym.Symbol(dico_Att[att][0])**int(i)
                expF = expF.replace(s, sym.Symbol(dico_Att[att][0]))
            dico_pw[dico_Att[att][0]] = 1
    return expF,dico_pw


"""
c1 = sym.Symbol('c1')
c4 = sym.Symbol('c4')
c5 = sym.Symbol('c5')
c6 = sym.Symbol('c6')
c7 = sym.Symbol('c7')
c9 = sym.Symbol('c9')
expression = 1 - c1*(1-c4)*(1-c5*(1-c6*(1-c7*(1-c9))))
res = sym.expand(expression)
print(res)


c2 = sym.Symbol('c2')
c3 = sym.Symbol('c3')
c8 = sym.Symbol('c8')
expression2 = 1 - c2*(1-c7*(1-c9))
res2 = sym.expand(expression2)
print(res2)

expression3 = 1 - c3*(1-c8*(1-c9))
res3 = sym.expand(expression3)
print(res3)

print('ici')
resfin = sym.expand(expression*expression2*expression3)
#resfin= 3*c1*c2*c3*c4*c5*c6*c7**2*c8*c9**2
print(resfin)

#sres = str(resfin)
#x = sres.replace("**3","")
#for i in range(2,4):
#s = sym.Symbol("**2")
#s2 = sym.Symbol("")
#resfin = resfin.replace(s, s2)
#s = sym.Symbol('c9')**int(3)
#resfin = resfin.replace(c9,sym.Symbol("XXX"))
#s = sym.Symbol("XXX")**int(3)
#resfin = resfin.replace(s,c9)

#resfin = resfin.replace(s, c7)
#print(resfin)

x = 0.1*0.2*0.3*0.4*0.5*0.6*0.7 *0.8*0.9  - 3*0.1*0.2*0.3*0.4*0.5*0.6*0.7 *0.8*0.9  + 3*0.1*0.2*0.3*0.4*0.5*0.6*0.7 *0.8*0.9 - 0.1*0.2*0.3*0.4*0.5*0.6*0.7 *0.8 + 0.1*0.2*0.3*0.4*0.5*0.6*0.7 *0.9  - 2*0.1*0.2*0.3*0.4*0.5*0.6*0.7 *0.9 + 0.1*0.2*0.3*0.4*0.5*0.6*0.7  + 2*0.1*0.2*0.3*0.4*0.5*0.6*0.7*0.8*0.9  - 4*0.1*0.2*0.3*0.4*0.5*0.6*0.7*0.8*0.9 + 2*0.1*0.2*0.3*0.4*0.5*0.6*0.7*0.8 + 2*0.1*0.2*0.3*0.4*0.5*0.6*0.7*0.9 - 2*0.1*0.2*0.3*0.4*0.5*0.6*0.7 + 0.1*0.2*0.3*0.4*0.5*0.6*0.8*0.9 - 0.1*0.2*0.3*0.4*0.5*0.6*0.8 + 0.1*0.2*0.3*0.4*0.5*0.6 - 0.1*0.2*0.3*0.4*0.5*0.7*0.8*0.9  + 2*0.1*0.2*0.3*0.4*0.5*0.7*0.8*0.9 - 0.1*0.2*0.3*0.4*0.5*0.7*0.8 - 0.1*0.2*0.3*0.4*0.5*0.7*0.9 + 0.1*0.2*0.3*0.4*0.5*0.7 - 0.1*0.2*0.3*0.4*0.5*0.8*0.9 + 0.1*0.2*0.3*0.4*0.5*0.8 - 0.1*0.2*0.3*0.4*0.5 + 0.1*0.2*0.3*0.4*0.7*0.8*0.9  - 2*0.1*0.2*0.3*0.4*0.7*0.8*0.9 + 0.1*0.2*0.3*0.4*0.7*0.8 + 0.1*0.2*0.3*0.4*0.7*0.9 - 0.1*0.2*0.3*0.4*0.7 + 0.1*0.2*0.3*0.4*0.8*0.9 - 0.1*0.2*0.3*0.4*0.8 + 0.1*0.2*0.3*0.4 - 0.1*0.2*0.3*0.5*0.6*0.7 *0.8*0.9  + 3*0.1*0.2*0.3*0.5*0.6*0.7 *0.8*0.9  - 3*0.1*0.2*0.3*0.5*0.6*0.7 *0.8*0.9 + 0.1*0.2*0.3*0.5*0.6*0.7 *0.8 - 0.1*0.2*0.3*0.5*0.6*0.7 *0.9  + 2*0.1*0.2*0.3*0.5*0.6*0.7 *0.9 - 0.1*0.2*0.3*0.5*0.6*0.7  - 2*0.1*0.2*0.3*0.5*0.6*0.7*0.8*0.9  + 4*0.1*0.2*0.3*0.5*0.6*0.7*0.8*0.9 - 2*0.1*0.2*0.3*0.5*0.6*0.7*0.8 - 2*0.1*0.2*0.3*0.5*0.6*0.7*0.9 + 2*0.1*0.2*0.3*0.5*0.6*0.7 - 0.1*0.2*0.3*0.5*0.6*0.8*0.9 + 0.1*0.2*0.3*0.5*0.6*0.8 - 0.1*0.2*0.3*0.5*0.6 + 0.1*0.2*0.3*0.5*0.7*0.8*0.9  - 2*0.1*0.2*0.3*0.5*0.7*0.8*0.9 + 0.1*0.2*0.3*0.5*0.7*0.8 + 0.1*0.2*0.3*0.5*0.7*0.9 - 0.1*0.2*0.3*0.5*0.7 + 0.1*0.2*0.3*0.5*0.8*0.9 - 0.1*0.2*0.3*0.5*0.8 + 0.1*0.2*0.3*0.5 - 0.1*0.2*0.3*0.7*0.8*0.9  + 2*0.1*0.2*0.3*0.7*0.8*0.9 - 0.1*0.2*0.3*0.7*0.8 - 0.1*0.2*0.3*0.7*0.9 + 0.1*0.2*0.3*0.7 - 0.1*0.2*0.3*0.8*0.9 + 0.1*0.2*0.3*0.8 - 0.1*0.2*0.3 - 0.1*0.2*0.4*0.5*0.6*0.7 *0.9  + 2*0.1*0.2*0.4*0.5*0.6*0.7 *0.9 - 0.1*0.2*0.4*0.5*0.6*0.7  - 2*0.1*0.2*0.4*0.5*0.6*0.7*0.9 + 2*0.1*0.2*0.4*0.5*0.6*0.7 - 0.1*0.2*0.4*0.5*0.6 + 0.1*0.2*0.4*0.5*0.7*0.9 - 0.1*0.2*0.4*0.5*0.7 + 0.1*0.2*0.4*0.5 - 0.1*0.2*0.4*0.7*0.9 + 0.1*0.2*0.4*0.7 - 0.1*0.2*0.4 + 0.1*0.2*0.5*0.6*0.7 *0.9  - 2*0.1*0.2*0.5*0.6*0.7 *0.9 + 0.1*0.2*0.5*0.6*0.7  + 2*0.1*0.2*0.5*0.6*0.7*0.9 - 2*0.1*0.2*0.5*0.6*0.7 + 0.1*0.2*0.5*0.6 - 0.1*0.2*0.5*0.7*0.9 + 0.1*0.2*0.5*0.7 - 0.1*0.2*0.5 + 0.1*0.2*0.7*0.9 - 0.1*0.2*0.7 + 0.1*0.2 - 0.1*0.3*0.4*0.5*0.6*0.7*0.8*0.9  + 2*0.1*0.3*0.4*0.5*0.6*0.7*0.8*0.9 - 0.1*0.3*0.4*0.5*0.6*0.7*0.8 - 0.1*0.3*0.4*0.5*0.6*0.7*0.9 + 0.1*0.3*0.4*0.5*0.6*0.7 - 0.1*0.3*0.4*0.5*0.6*0.8*0.9 + 0.1*0.3*0.4*0.5*0.6*0.8 - 0.1*0.3*0.4*0.5*0.6 + 0.1*0.3*0.4*0.5*0.8*0.9 - 0.1*0.3*0.4*0.5*0.8 + 0.1*0.3*0.4*0.5 - 0.1*0.3*0.4*0.8*0.9 + 0.1*0.3*0.4*0.8 - 0.1*0.3*0.4 + 0.1*0.3*0.5*0.6*0.7*0.8*0.9  - 2*0.1*0.3*0.5*0.6*0.7*0.8*0.9 + 0.1*0.3*0.5*0.6*0.7*0.8 + 0.1*0.3*0.5*0.6*0.7*0.9 - 0.1*0.3*0.5*0.6*0.7 + 0.1*0.3*0.5*0.6*0.8*0.9 - 0.1*0.3*0.5*0.6*0.8 + 0.1*0.3*0.5*0.6 - 0.1*0.3*0.5*0.8*0.9 + 0.1*0.3*0.5*0.8 - 0.1*0.3*0.5 + 0.1*0.3*0.8*0.9 - 0.1*0.3*0.8 + 0.1*0.3 + 0.1*0.4*0.5*0.6*0.7*0.9 - 0.1*0.4*0.5*0.6*0.7 + 0.1*0.4*0.5*0.6 - 0.1*0.4*0.5 + 0.1*0.4 - 0.1*0.5*0.6*0.7*0.9 + 0.1*0.5*0.6*0.7 - 0.1*0.5*0.6 + 0.1*0.5 - 0.1 + 0.2*0.3*0.7*0.8*0.9  - 2*0.2*0.3*0.7*0.8*0.9 + 0.2*0.3*0.7*0.8 + 0.2*0.3*0.7*0.9 - 0.2*0.3*0.7 + 0.2*0.3*0.8*0.9 - 0.2*0.3*0.8 + 0.2*0.3 - 0.2*0.7*0.9 + 0.2*0.7 - 0.2 - 0.3*0.8*0.9 + 0.3*0.8 - 0.3 + 1
print(x)
#Pow(x, 2)

#c1*c2*c3*c4*c5*c6*c7**2*c8*c9**3 - c1*c2*c3*c4*c5*c6*c7**2*c8*c9**2 + c1*c2*c3*c5*c6*c7**2*c8*c9**3 - c1*c2*c3*c5*c6*c7**2*c8*c9**2
#c1*c2*c3*c4*c5*c6*(c7-c7**2)*c8*(c9-c9**3) - c1*c2*c3*c4*c5*c6*(c7-c7**2)*c8*(c9-c9**2) + c1*c2*c3*c5*c6*(c7-c7**2)*c8*(c9-c9**3) - c1*c2*c3*c5*c6*(c7-c7**2)*c8*(c9-c9**2)
#e = (c1*c4*c5*c6*c7*c9 - c1*c4*c5*c6*c7 - c1*c5*c6*c7 + c1*c5*c6*c7*c9) * (-c2*c7*c9) * (-c3*c8*c9)
#re = sym.expand(e)
#print(re)


#x = sym.Symbol('x')
#y = sym.Symbol('y')
#z = sym.Symbol('z')

#form = 1 -w -x -y -z + x*y + x*z + y*z +x*w + y*w + z*w - x*y*z - x*y*w - y*z*w + w*x*y*z
#print(sym.factor(form))

#exfull = expression * expression2 * expression3
#res0 = sym.expand(exfull)
#print(res0)
#0.1*0.2*0.3*0.4*0.5*0.6*(0.7-0.7**2)*0.8*(0.9-0.9**3) - 0.1*0.2*0.3*0.4*0.5*0.6*(0.7-0.7**2)*0.8*(0.9-0.9**2) = 0.00000979776
# + 0.1*0.2*0.3*0.5*0.6*(0.7-0.7**2)*0.8*(0.9-0.9**3) - 0.1*0.2*0.3*0.5*0.6*(0.7-0.7**2)*0.8*(0.9-0.9**2) = 0.0000244944
"""

############################################################


def Fast(goal, dico_Arg, dico_Att, dico_lvl,dico_pw):
    ldico_att_in = list_dico_Att_in(dico_Arg, dico_Att)
    dico_att_in = dlist_Att_to_Arg(ldico_att_in)
    dico_lvl = level_Arg(goal, 1, dico_att_in, dico_lvl)
    liste_lvl = order_level(dico_lvl)
    print(liste_lvl)
    liste_lvl = liste_lvl[1:]
    print(liste_lvl)

    formula = 1
    for att in ldico_att_in[goal]:
        formula = formula * (1- sym.Symbol(dico_Att[att][0]) * float(-dico_Att[att][2]))
        dico_pw[dico_Att[att][0]] += 1

    print(f'List of Levels: {liste_lvl}')
    for arg in liste_lvl:
        formula, dico_pw = develop(formula, arg[0], ldico_att_in, dico_Att, dico_pw)

    return formula

start = time.time()
proba = Fast("a", dico_Arg, dico_Att, dico_lvl, dico_pw)
end = time.time()
elapsed = end - start

print(f'Probability of a = {proba}')
print(f'Time: {elapsed:.5}s\n')