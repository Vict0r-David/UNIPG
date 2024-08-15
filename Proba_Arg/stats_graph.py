#import math
import numpy  

import time

#from sage.all import *
#from sage.arith.power import generic_power

import matplotlib.pyplot as plt

import statistics

#import pickle



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

def dlist_Att_out_to_Arg(dico_list_att):
    dico = {}
    for arg,l in dico_list_att.items():
        dico[arg] = []
        for att in l:
            l2 = att.split("->")
            dico[arg].append(l2[1])
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

############################################### Acyclic Case ##############################################

# Preprocessing

def level_Arg(start, lvl, ldico_att_in,dico_lvl):
    for att in ldico_att_in[start]:
        if dico_lvl[att] < lvl:
            dico_lvl[att] = lvl
        level_Arg(att, lvl+1, ldico_att_in,dico_lvl)
    return dico_lvl

#def init_dico_lvl(dico_Arg,dico_lvl):
#    for arg in dico_Arg:
#        dico_lvl[arg] = 0
#    return dico_lvl


def order_level(start,dico_lvl):
    liste = []
    add = False
    for arg,val in dico_lvl.items():
            if val > 0 or arg == start:
                liste.append([arg,val])
    heapSort(liste)
    #size = len(liste)
    #quickSort(liste, 0, size - 1)
    return liste


############################################################    



#############################################################################################
################################ INUTILE dans cette version #################################

def list_dico_Att_out(dico_Arg, dico_Att):
    ldico_att_out = {}
    for id_arg in dico_Arg:
        ldico_att_out[id_arg] = []
    for id,att in dico_Att.items():
        ldico_att_out[att[0]] += [id]
    return ldico_att_out


def get_paths_back2(node, ldico_att_in, dico_att, paths=None, current_path=None, back_arg=None):
    if paths is None:
        paths = []
    if current_path is None:
        current_path = []
    if back_arg is None:
        back_arg = set()
    else:
        back_arg.add(node)

    current_path.append(node)
    if ldico_att_in[node] == []: 
        paths.append(current_path)
    else:
        children = []
        for att in ldico_att_in[node]:
            children.append(dico_att[att][0])
        for child in children:
            get_paths_back2(child, ldico_att_in, dico_att, paths, list(current_path),back_arg)
    return paths,back_arg


def build_dico_paths2(goal, ldico_att_in, dico_att):
    list_paths,back_att = get_paths_back2(goal, ldico_att_in, dico_att, paths=None, current_path=None,back_arg=None)
    #print(f"list paths: {list_paths}")
    #print(f"nb paths: {len(list_paths)}")
    back_att.add(goal)
    dico_paths = {}
    for path in list_paths:
        for i in range(0,len(path)):
            if str(path[i]) not in dico_paths:
                dico_paths[path[i]] = set(path[i:])
            else:
                dico_paths[str(path[i])].update(set(path[i:]))
    return dico_paths,back_att

def dependant_arg2(dico_att_out,back_att):
    dep = set()
    for a in back_att:
        cpt = 0
        for att in dico_att_out[a]:
            if att in back_att:
                cpt +=1
            if cpt > 1:
                dep.add(a)
                break
    return dep


def init_dico_lvl(dico_Arg,dico_lvl):
    for arg in dico_Arg:
        dico_lvl[arg] = 0
    return dico_lvl
    

#### GESTION DU GRAPHE  - INITIALISATION ####

#output = "./50_75_1500/output_MCN_DAG_50_75_1500.txt"
#f_out = open(output,"w")

liste_in = []
liste_out = []

"""
pickle.dump([],open("i_max_50_75_1000", "wb"))
pickle.dump([],open("i_sum_50_75_1000", "wb"))
pickle.dump([],open("i_time_50_75_1000", "wb"))

moyenne_dep = 3.9575152880592213
ecart_pop_dep = 3.6588681840377144

moyenne_kmin = 3.3109800319102645
ecart_pop_kmin = 1.5407842709307518

moyenne_kmax = 6.764831020644974
ecart_pop_kmax = 1.9790047503690167

"""

#f_out.write("Graph ")
#f_out.write(str(i))

nb_graph = 0

list_dep = []
list_kmin = []
list_kmax = []

dico_dep = {}
dico_kmin = {}
dico_kmax = {}

for i in range(1,1501):
    print("=================================")
    print("=================================")
    print(i)
    print("=================================")
    print("=================================")
    file_AF = "./50_75_1500/DAG_50_75_"+str(i)+".txt"
    #file_AF = ".\SCN_1000\SCN_1000_"+str(i)+".txt"
    AF = open(file_AF,"r")
    dico_Arg = {}
    dico_Att = {}
    dico_lvl = {}
    dico_pw = {}

    slist_dep = []
    slist_kmin = []
    slist_kmax = []

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

    

    dico_list_att_in = list_dico_Att_in(dico_Arg, dico_Att)
    d_att_in = dlist_Att_to_Arg(dico_list_att_in)

    dico_list_att_out = list_dico_Att_out(dico_Arg, dico_Att)
    d_att_out = dlist_Att_to_Arg(dico_list_att_out)

    avg_dep = 0
    avg_kmin = 0
    avg_kmax = 0

    nb_arg = 0

    for a in dico_Arg:
        print(a)
        if len(d_att_in[a]) >= 0:
            liste_in.append(len(d_att_in[a]))
        if len(d_att_out[a]) >= 0:
            liste_out.append(len(d_att_out[a]))

        d_lvl = init_dico_lvl(dico_Arg,{})
        dico_lvl = level_Arg(a, 1, d_att_in, d_lvl)    
        liste_lvl = order_level(a,dico_lvl)
        dico_deg = {}
        i = len(liste_lvl)-1

        dico_paths2, back_att = build_dico_paths2(a, dico_list_att_in, dico_Att)
        merge_dep = dependant_arg2(d_att_out, back_att)   

        avg_dep += len(merge_dep)
        list_dep.append(len(merge_dep))

        if len(merge_dep) > 0:
            mini = 100
            maxi = 0
            for d in merge_dep:
                dep_att = 0
                for att in d_att_out[d]:
                    if att in back_att:
                        dep_att += 1
                if dep_att > maxi:
                    maxi = dep_att
                if dep_att < mini:
                    mini = dep_att
            #if mini == 100:
            #    mini = 2
            #if maxi == 0:
            #    maxi = 2

            avg_kmin += mini
            avg_kmax += maxi

            list_kmin.append(mini)
            list_kmax.append(maxi)

            
            slist_kmin.append(mini)
            slist_kmax.append(maxi)
        slist_dep.append(len(merge_dep))

    dico_dep[i] = slist_dep
    dico_kmin[i] = slist_kmin
    dico_kmax[i] = slist_kmax


moy_dep = 0
moy_ecart = 0
nb = 0
for key,liste in dico_dep.items():
    moy_dep += statistics.mean(liste)
    moy_ecart += statistics.pstdev(liste)
    nb += 1

print(f"moyenne_dep par graphe = {moy_dep/nb}")
print(f"ecart_pop_dep par graphe = {moy_ecart/nb}\n")

######

moy_kmin = 0
moy_ecart = 0
nb = 0
for key,liste in dico_kmin.items():
    moy_kmin += statistics.mean(liste)
    moy_ecart += statistics.pstdev(liste)
    nb += 1

print(f"moyenne_kmin par graphe = {moy_kmin/nb}")
print(f"ecart_pop_kmin par graphe = {moy_ecart/nb}\n")

######

moy_kmax = 0
moy_ecart = 0
nb = 0
for key,liste in dico_kmax.items():
    moy_kmax += statistics.mean(liste)
    moy_ecart += statistics.pstdev(liste)
    nb += 1

print(f"moyenne_kmax par graphe = {moy_kmax/nb}")
print(f"ecart_pop_kmax par graphe = {moy_ecart/nb}\n")


#print(f"list dep = {list_dep}")
#print(f"list kmin = {list_kmin}")
#print(f"list kmax = {list_kmax}")

moyenne_dep = statistics.mean(list_dep)
ecart_pop_dep = statistics.pstdev(list_dep)
print(f"moyenne_dep = {moyenne_dep}")
print(f"ecart_pop_dep = {ecart_pop_dep}\n")

moyenne_kmin = statistics.mean(list_kmin)
ecart_pop_kmin = statistics.pstdev(list_kmin)
print(f"moyenne_kmin = {moyenne_kmin}")
print(f"ecart_pop_kmin = {ecart_pop_kmin}\n")

moyenne_kmax = statistics.mean(list_kmax)
ecart_pop_kmax = statistics.pstdev(list_kmax)
print(f"moyenne_kmax = {moyenne_kmax}")
print(f"ecart_pop_kmax = {ecart_pop_kmax}\n")

dico = {}   
for elem in liste_in:
    if elem in dico:
        dico[elem] += 1
    else:
        dico[elem] = 1

proportion = []
for key,val in dico.items():
    proportion.append([key,val/len(liste_in)])
    
moyenne = statistics.mean(liste_in)
ecart_pop = statistics.pstdev(liste_in)

print(f"mean in = {moyenne}")
print(f"ecart pop in = {ecart_pop}")
print(f"dico in = {dico}")
print(f"proportion in = {proportion}")

"""
mean = 2.062290333816062
ecart pop = 1.1666852485282153
dico = {1: 21948, 2: 16735, 3: 9404, 4: 4293, 5: 1531, 7: 130, 6: 477, 9: 7, 8: 26}
proportion = [[1, 0.40], [2, 0.306], [3, 0.172], [4, 0.0786], [5, 0.028], [6, 0.0087]
proportion = [[1, 0.39], [2, 0.31], [3, 0.17], [4, 0.08], [5, 0.03], [6, 0.01],
"""

        
dico_out = {}   
for elem in liste_out:
    if elem in dico_out:
        dico_out[elem] += 1
    else:
        dico_out[elem] = 1

proportion_out = []
for key,val in dico_out.items():
    proportion_out.append([key,val/len(liste_out)])
    
moyenne_out = statistics.mean(liste_out)
ecart_pop_out = statistics.pstdev(liste_out)

print(f"mean out = {moyenne_out}")
print(f"ecart pop out = {ecart_pop_out}")
print(f"dico out = {dico_out}")
print(f"proportion out = {proportion_out}")

"""
mean out = 3.240860772621208
ecart pop out = 1.9617313855620024
dico out = {3: 6484, 5: 3685, 2: 7126, 4: 5124, 1: 7631, 8: 625, 6: 2242, 7: 1345, 9: 255, 10: 120, 11: 49, 12: 16, 13: 5, 16: 1, 14: 5}
proportion out = [[3, 0.18678881110823034], [5, 0.10615619508541468], [2, 0.20528332325065538], [4, 0.14761040532365396], [1, 0.21983118716331057], [8, 0.01800478207011782], [6, 0.06458675424192666], [7, 0.038746291014893554], [9, 0.007345951084608072], [10, 0.003456918157462622], [11, 0.0014115749142972374], [12, 0.0004609224209950163], [13, 0.00014403825656094258], [16, 2.880765131218852e-05], [14, 0.00014403825656094258]]
proportion out = [[1, 0.220], [2, 0.205], [3, 0.187], [4, 0.148], [5, 0.106], [6, 0.065], [7, 0.039], [8, 0.018], [9, 0.007] ]
proportion out = [[1, 0.22], [2, 0.21], [3, 0.19], [4, 0.15], [5, 0.11], [6, 0.07], [7, 0.04], [8, 0.02], [9, 0.01] ]
"""


"""
Statistique Attaque entrante pour 1 arg:
mean in = 2.06
ecart pop in = 1.17
proportion in = [[1, 0.39], [2, 0.31], [3, 0.17], [4, 0.08], [5, 0.03], [6, 0.01]]

Statistique Attaque sortante pour 1 arg:
mean out = 3.24
ecart pop out = 1.96
proportion out = [[1, 0.22], [2, 0.20], [3, 0.19], [4, 0.15], [5, 0.10], [6, 0.06], [7, 0.04], [8, 0.02], [9, 0.01]]
"""
