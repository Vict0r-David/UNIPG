# Victor DAVID
# Created 03/10/2022

#import math
import numpy  

import time

from sage.all import *
from sage.arith.power import generic_power

import matplotlib.pyplot as plt

import pickle

import Approx as MC


#file_AF = "./50_75_1000/DAG_50_75_2.txt"

"""

#### GESTION DU GRAPHE  - INITIALISATION ####

#file_AF = ".\Old_test\AF5_3.txt"
#file_AF = "./DATA/DAG_50_80/DAG_50_80_10.txt"
#file_AF = "./DATA/DAG_50_90/DAG_50_90_10.txt"
#file_AF = "./DATA/DAG_50_100/DAG_50_100_1.txt"

#file_AF = "./DAG-test.txt"
#file_AF = "./small-test.txt"
#file_AF = "./test-paper.txt"
#file_AF = "./diamond.txt"

#file_AF = "./DAG-Test-small.txt"      
#file_AF = "./diamond1.txt"     

#file_AF = "./test1.txt"
#file_AF = "./test2.txt"
#file_AF = "./test3.txt"
#file_AF = "./verif.txt"

#file_AF = ".\SCN_500\SCN_500_1.txt"
#file_AF = ".\Old_test\AF5_33.txt"
AF = open(file_AF,"r")
dico_Arg = {}
dico_Att = {}
dico_lvl = {}
dico_pw = {}
dico_pw_att = {}

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
        dico_pw_att[id] = 0


#print(dico_Arg)
#print(dico_Att)

#AF_1.txt
#dico_Arg = {"a":1, "b":1, "c":1, "d":1}
#dico_Att = {"a->b":["a","b",-0.6], "b->c": ["b","c",-0.3], "d->b":["d","b",-0.2], "a->d":["a","d",-0.5], "c->a":["c","a",-0.2]}

#AF_2.txt
#dico_Arg = {"a":1, "b":1, "c":1, "d":1, "e":1}
#dico_Att = {"a->c": ["a","c",-0.3], "b->c": ["b","c",-0.9], "c->e": ["c","e",-0.4], "d->e": ["d","e",-0.3]}

"""

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



#############################################################################################
################################ INUTILE dans cette version #################################

def list_dico_Att_out(dico_Arg, dico_Att):
    ldico_att_out = {}
    for id_arg in dico_Arg:
        ldico_att_out[id_arg] = []
    for id,att in dico_Att.items():
        ldico_att_out[att[0]] += [id]
    return ldico_att_out

"""
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
"""

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
"""
def build_dico_paths(goal, ldico_att_in, dico_att):
    list_paths = get_paths_back(goal, ldico_att_in, dico_att, paths=None, current_path=None)
    #print(f"list paths: {list_paths}")
    #print(f"nb paths: {len(list_paths)}")
    dico_paths = {}
    for path in list_paths:
        for i in range(0,len(path)):
            if str(path[i]) not in dico_paths:
                dico_paths[path[i]] = set(path[i:])
            else:
                dico_paths[str(path[i])].update(set(path[i:]))
    return dico_paths


def dependant_arg(goal,dico_att_in,dico_att_out,dico_paths):
    back_att = set()
    for arg,path in dico_paths.items():
        back_att.update(path)
    dep = set()
    for a in back_att:
        cpt = 0
        for b,path in dico_paths.items():
            if a in path:
                cpt += 1
            if cpt > 1:
                cpt2 = 0
                for c in dico_att_out[a]:
                    if c in back_att:
                        cpt2 += 1
                    if cpt2 > 1:
                        dep.add(a)
                        break
                break
    return dep
"""

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

############## Function to understand complexity ############################# 

"""
def conjunction_with_larger(goal, dico_att_in, dico_deg, dep):
    dico_conj = {}
    moyenne = 0
    maxi = 0
    list_arg = [goal]
    dejavu = [goal]

    while len(list_arg) > 0:
        symb_att = 0
        for att in dico_att_in[list_arg[0]]:
            if att not in dejavu:
                list_arg.append(att)
                dejavu.append(att)
            if isinstance(dico_deg[att], sym.Basic) or att in dep:
                symb_att += 1
        if symb_att > 1:
            dico_conj[list_arg[0]] = symb_att
            moyenne += symb_att
            if symb_att > maxi:
                maxi = symb_att
        del list_arg[0]

    nb = len(dico_conj)
    moyenne = moyenne/nb

    print(f"nb conj = {nb}")
    print(f"max larger = {maxi}")
    print(f"avg larger = {moyenne}")

    return dico_conj
"""

def nbTermes(liste_lvl, dico_att_in, dico_deg, dep):
    dico_term = {}
    dico_term_symb = {}
    i = len(liste_lvl)-1

    list_term = []
    #som_terms = 0

    #Termes Symboliques = Prod_VarDist 2 * [0.5 + sommeNBVar 0.5]

    while i >= 0: 
        l_att = dico_att_in[str(liste_lvl[i][0])]
        dico_term[str(liste_lvl[i][0])] = 1
        dico_term_symb[str(liste_lvl[i][0])] = {}
        if type(dico_deg[str(liste_lvl[i][0])]) == sage.symbolic.expression.Expression:
            for att in l_att:
                if att in dep:
                    if str(att) not in dico_term_symb[str(liste_lvl[i][0])].keys():
                        dico_term_symb[str(liste_lvl[i][0])][str(att)] = 1
                    else :
                        dico_term_symb[str(liste_lvl[i][0])][str(att)] += 1

                elif type(dico_deg[att]) == sage.symbolic.expression.Expression:
                    for key,val in dico_term_symb[str(att)].items():
                        if str(key) in dico_term_symb[str(liste_lvl[i][0])]:
                            dico_term_symb[str(liste_lvl[i][0])][str(key)] += val
                        else:
                            dico_term_symb[str(liste_lvl[i][0])][str(key)] = val

            for key,val in dico_term_symb[str(liste_lvl[i][0])].items():
                dico_term[str(liste_lvl[i][0])] = dico_term[str(liste_lvl[i][0])] * 2 * (0.5 + 0.5 * val)
        i-=1

    #print(f"dico term {dico_term}")

    if len(dep) > 0:
        j = 0
        terms =  dico_term[str(liste_lvl[0][0])]
        #som_terms = terms
        list_term.append(terms)
        a = str(liste_lvl[0][0])

        # keep only the dependant argument and thanks to list_lvl they are ordered increasingly
        while j < len(liste_lvl):
            if liste_lvl[j][0] not in dep:
                del liste_lvl[j]
                j -=1
            j +=1

        #print(f"dts: {dico_term_symb}")

        for arg,lvl in liste_lvl:   

            for key,val in dico_term_symb[arg].items():
                if key in dico_term_symb[a]:
                    dico_term_symb[a][key] += val
                else:
                    dico_term_symb[a][key] = val

            del dico_term_symb[a][arg]

            t = 1
            for key,val in dico_term_symb[a].items():
                t = t * 2 * (0.5 + 0.5 * val)
            #som_terms = som_terms + t

            list_term.append(t)
        
        #print(f"list terms = {list_term}")
        #print(f"max terms = {max(list_term)}")
        #print(f"sum terms = {som_terms}\n")

    if len(list_term) == 0:
        list_term.append(0)

    return dico_term,list_term,max(list_term)

############## ############## ############## ############## ############## 

"""
def fastMCN(goal,dico_Arg,dico_Att,d_lvl):
    ldico_att_in = list_dico_Att_in(dico_Arg, dico_Att)
    dico_att_in = dlist_Att_to_Arg(ldico_att_in)
    dico_att_out = list_dico_Att_out(dico_Arg, dico_Att)
    d_att_out =  dlist_Att_out_to_Arg(dico_att_out)

    start = time.time()

    dico_lvl = level_Arg(goal, 1, dico_att_in, d_lvl)
    liste_lvl = order_level(goal,dico_lvl)
    dico_deg = {}
    i = len(liste_lvl)-1
    #dico_paths = build_dico_paths(goal, ldico_att_in, dico_Att)
    dico_paths,back_att = build_dico_paths2(goal, ldico_att_in, dico_Att)
    
    #merge_dep = dependant_arg(goal, dico_att_in,d_att_out, dico_paths)
    merge_dep = dependant_arg2(dico_att_out,back_att)    
    
    symbo = False
    while i >= 0: 
        l_att = ldico_att_in[str(liste_lvl[i][0])]
        deg = 1
        for att in l_att:
                if dico_Att[att][0] in merge_dep :
                    #deg = deg * (1-(sym.Symbol(dico_Att[att][0])* -dico_Att[att][2]))
                    deg = deg * (1-(var(dico_Att[att][0])* -dico_Att[att][2]))
                    symbo = True
                else:
                    deg = deg * (1-(dico_deg[dico_Att[att][0]]*- dico_Att[att][2])) 
        
        dico_deg[str(liste_lvl[i][0])] = deg
        
        i-=1

    liste = list(liste_lvl)

    if len(merge_dep) > 0:
            
            j = 0
            # keep only the dependant argument and thanks to list_lvl they are ordered increasingly
            while j < len(liste_lvl):
                if liste_lvl[j][0] not in merge_dep:
                    del liste_lvl[j]
                    j -=1
                j +=1

            for arg,lvl in liste_lvl:   

                #Expand
                #start = time.time()
                dico_deg[goal] = dico_deg[goal].expand()
                #end = time.time()
                #elapsed = end - start
                #print(f"{arg} exp : {elapsed}")

                #Delete
                for argu,lvl in liste_lvl:
                    for i in range(2,dico_deg[goal].degree(var(str(argu)))+1):
                        dico_deg[goal] = dico_deg[goal].subs(generic_power(var(str(argu)), i)==var(str(argu)))

                #Replace 
                dico_deg[goal] = dico_deg[goal].subs(var(str(arg))==dico_deg[arg]) 

    end = time.time()
    elapsed = end - start
    #print(f"Time : {elapsed}")

    dico_term, list_term, max_list_term = nbTermes(liste, dico_att_in, dico_deg, merge_dep)

    return dico_deg[goal], elapsed, list_term, max_list_term
"""

def fastMCN(goal,dico_Arg,dico_Att,d_lvl):
    ldico_att_in = list_dico_Att_in(dico_Arg, dico_Att)
    dico_att_in = dlist_Att_to_Arg(ldico_att_in)
    dico_att_out = list_dico_Att_out(dico_Arg, dico_Att)
    d_att_out =  dlist_Att_out_to_Arg(dico_att_out)

    start = time.time()

    dico_lvl = level_Arg(goal, 1, dico_att_in, d_lvl)
    liste_lvl = order_level(goal,dico_lvl)
    dico_deg = {}
    i = len(liste_lvl)-1

    #####
    dico_paths2, back_att = build_dico_paths2(goal, ldico_att_in, dico_Att)
    
    merge_dep = dependant_arg2(d_att_out, back_att)
    
    symbo = False
    while i >= 0: 
        l_att = ldico_att_in[str(liste_lvl[i][0])]
        deg = 1
        for att in l_att:
                if dico_Att[att][0] in merge_dep :
                    #deg = deg * (1-(sym.Symbol(dico_Att[att][0])* -dico_Att[att][2]))
                    deg = deg * (1-(var(dico_Att[att][0])* -dico_Att[att][2]))
                    symbo = True
                else:
                    deg = deg * (1-(dico_deg[dico_Att[att][0]]*- dico_Att[att][2])) 
        dico_deg[str(liste_lvl[i][0])] = deg
        
        #if type(deg) == sage.symbolic.expression.Expression:
        #    dico_deg[str(liste_lvl[i][0])] = deg.expand()
        #else:
        #    dico_deg[str(liste_lvl[i][0])] = deg

        i-=1

    liste = list(liste_lvl)
    if len(merge_dep) > 0:
            
            j = 0
            # keep only the dependant argument and thanks to list_lvl they are ordered increasingly
            while j < len(liste_lvl):
                if liste_lvl[j][0] not in merge_dep:
                    del liste_lvl[j]
                    j -=1
                j +=1

            for arg,lvl in liste_lvl:   

                #Expand
                dico_deg[goal] = dico_deg[goal].expand()

                #Delete
                for argu,lvl in liste_lvl:

                    for i in range(2,dico_deg[goal].degree(var(str(argu)))+1):
                        dico_deg[goal] = dico_deg[goal].subs(generic_power(var(str(argu)), i)==var(str(argu)))

                #Replace 
                dico_deg[goal] = dico_deg[goal].subs(var(str(arg))==dico_deg[arg]) 

    end = time.time()
    elapsed = end - start
    #print(f"Time : {elapsed}")

    dico_term, list_term, max_list_term = nbTermes(liste, dico_att_in, dico_deg, merge_dep)

    #print(f"dico term {dico_term}")
    #print(f"list term {list_term}")
    #print(f"max_list_term {max_list_term}")
    #print(f"som_terms {som_terms}")

    #conjunction_with_larger(goal, dico_att_in, dico_deg, merge_dep)

    return dico_deg[goal], elapsed, list_term, max_list_term

############################################################

"""
start = time.time()
p = fastMCN("a",dico_Arg,dico_Att,dico_lvl)
end = time.time()
elapsed = end - start

print(f"Probability of a (ALL MCN algo) = {p}")
print(f"Time (ALL MCN algo) = {elapsed}\n")
"""



###########################################################################################################
############################################# Experimentation #############################################

def init_dico_lvl(dico_Arg,dico_lvl):
    for arg in dico_Arg:
        dico_lvl[arg] = 0
    return dico_lvl

"""
start = time.time()
p = fastMCN("a16",dico_Arg,dico_Att,dico_lvl)
end = time.time()
elapsed = end - start

print(f"Probability of a (ALL MCN algo) = {p}")
print(f"Time (ALL MCN algo) = {elapsed}\n")

"""

#### GESTION DU GRAPHE  - INITIALISATION ####

output = "./50_75_1000/RAPIDE_output_MCN_DAG_50_75_1000.txt"
#output = ".\SCN_1000\output_M_SCN_1000.txt"
f_out = open(output,"w")
avg_total_time = 0
avg_total_nodes = 0
avg_total_edges = 0
total_max = 0

tab_Max = []
#tab_Sum = []
tab_Time = []

dico_FAST_MCN = {}
dico_MC_MCN = {}
#dico_diff = {}

List_diff = []
List_error = []

avg_diff = 0
avg_error = 0

cpt = 0

compteur = 0

cpt2 = 0
cptProb_0 = 0
cpt_0 = 0

f_out.write("Graph ")
f_out.write(str(i))
#1501
for i in range(1,1001):
    print("=================================")
    print("=================================")
    print(i)
    print("=================================")
    print("=================================")
    file_AF = "./50_75_1000/DAG_50_75_"+str(i)+".txt"
    #file_AF = ".\SCN_1000\SCN_1000_"+str(i)+".txt"
    AF = open(file_AF,"r")
    dico_Arg = {}
    dico_Att = {}
    dico_lvl = {}
    dico_pw = {}

    dico_FAST_MCN[i] = {} 
    dico_MC_MCN[i] = {} 
    #dico_diff[i] = {}

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

    max_time = 0
    min_time = 999999999999999999
    liste_output = []


    for a in dico_Arg:
        print(a)
        #start = time.time()
        proba, time_MCN, list_term, max_list_term = fastMCN(str(a),dico_Arg,dico_Att,dico_lvl)
        #end = time.time()
        #elapsed = end - start
        tab_Max.append(max_list_term)
        #tab_Sum.append(som_terms)
        tab_Time.append(round(time_MCN,4))

        if time_MCN < 0.00001:
            time_MCN = 0.00001


        dico_FAST_MCN[i][a] = [proba,time_MCN] 

        
        approximation = MC.approx(str(a), time_MCN, dico_Att, dico_Arg)
        dico_MC_MCN[i][a] = [approximation,time_MCN]

        if proba == approximation:
            cpt2+=1

        diff = math.sqrt((proba - approximation)**2)

        #List_diff.append(diff)
        #avg_diff += diff

        if proba == 0 and approximation != 0:
            #List_error.append(100)
            #avg_error += 100
            cptProb_0 +=1

        elif proba == 0 and approximation == 0:
            List_error.append(0)
            avg_error += 0
            List_diff.append(0)
            avg_diff += 0
            cpt_0 += 1

        else:
            #if float(diff/proba) == 1:
            if approximation == 0:
                compteur += 1

            else: 
                List_error.append(float(diff/proba)*100)
                avg_error += float(diff/proba)*100
                List_diff.append(diff)
                avg_diff += diff
        
        cpt += 1

        liste_output.append([proba,round(time_MCN,4),max_list_term])
        dico_lvl = init_dico_lvl(dico_Arg,dico_lvl)

        if time_MCN > max_time:
            max_time = time_MCN
        if time_MCN < min_time:
            min_time = time_MCN
        
        #print(proba,round(elapsed,4))

    #print(f'Probability of a = {proba}')
    #print(f'Time: {elapsed:.5}s\n')

        f_out.write("Prob arg(")
        f_out.write(str(a))
        f_out.write(") = ")
        f_out.write(str(proba))
        f_out.write(" and Time = ")
        f_out.write(str(round(time_MCN,4)))
        f_out.write(" max terms = ")
        f_out.write(str(max_list_term))
        #f_out.write(" sum terms = ")
        #f_out.write(str(som_terms))
        f_out.write("\n")

    f_out.write("Average_Time = ")
    avg = 0
    for pair in liste_output:
        avg += pair[1]
    avg = avg/len(liste_output)
    f_out.write(str(avg))
    nb_nodes = len(dico_Arg)
    nb_edges = len(dico_Att)
    f_out.write(" nb_nodes = ")
    f_out.write(str(nb_nodes))
    f_out.write(" nb_edges = ")
    f_out.write(str(nb_edges))
    f_out.write("\nMin_Time = ")
    f_out.write(str(min_time))
    f_out.write(" Max_Time = ")
    f_out.write(str(max_time))
    f_out.write("\n \n")

    if max_time > total_max:
        total_max = max_time

    avg_total_time += avg
    avg_total_nodes += nb_nodes
    avg_total_edges += nb_edges

    AF.close()

avg_total_time = avg_total_time/1000
avg_total_nodes = avg_total_nodes/1000
avg_total_edges = avg_total_edges/1000
f_out.write("Average Total Time = ")
f_out.write(str(avg_total_time))
f_out.write("\nAverage Total Nodes = ")
f_out.write(str(avg_total_nodes))
f_out.write("\nAverage Total Edges = ")
f_out.write(str(avg_total_edges))
f_out.write("\nMax Total Time for one argument = ")
f_out.write(str(total_max))
f_out.close()

#print(tab_Max)
#print(tab_Sum)
#print(tab_Time)


pickle.dump(tab_Max, open("max_50_75_1000_rapide", "wb"))
#pickle.dump(tab_Sum, open("sum_50_75_1500", "wb"))
pickle.dump(tab_Time, open("time_50_75_1000_rapide", "wb"))


avg_diff = avg_diff/len(List_diff)    
avg_error = avg_error/len(List_error)       

print(f"avg diff = {avg_diff}")
print(f"avg error = {avg_error}")
print(f"approximation == 0, = {compteur}")
print(f"proba == approximation, = {cpt2}")
print(f"cptProb_0 et approx not 0, = {cptProb_0}")
print(f"nb valeur error percentage = {len(List_error)}")



l_error_small = []
for err in List_error:
    if err <= 100 and err>=0:
        l_error_small.append(err)

#pickle.dump(dico_diff, open("Diff_SCN", "wb"))

dico_FAST_MCN = {}
dico_MC_MCN = {}

pickle.dump(List_diff, open("diff_50_75_1000_rapide", "wb"))
pickle.dump(List_error, open("error_50_75_1000_rapide", "wb"))
pickle.dump(l_error_small, open("error_borne_50_75_1000_rapide", "wb"))

pickle.dump(dico_FAST_MCN, open("DicoFAST_50_75_1000_rapide", "wb"))
pickle.dump(dico_MC_MCN, open("DicoMC_50_75_1000_rapide", "wb"))



figure = plt.figure(figsize = (10, 10))
plt.gcf().subplots_adjust(left = 0.2, bottom = 0.2,
                       right = 0.7, top = 0.7, 
                       wspace = 0.5, hspace = 0.7)

axes = figure.add_subplot(1, 2, 1)
plt.hist(List_diff,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
#plt.title("Histogram of 1000 MCN (50 arg, 75 att)")
plt.xlabel('Error difference')

#plt.yticks([0,1000,3000,5000,10000,15000,20000,25000,30000])
plt.ylabel('Number of occurrence')


axes = figure.add_subplot(1, 2, 2)
#plt.hist(List_error,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
plt.hist(l_error_small,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
#plt.title("Histogram of 1000 MCN (50 arg, 75 att)")

#plt.yticks([0,1000,3000,5000,10000,15000,20000])
plt.xlabel('Error percentage')
plt.ylabel('Number of occurrence')

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

