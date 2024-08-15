# Victor DAVID
# Created 03/10/2022

#import math
import numpy  

import time

from sage.all import *
from sage.arith.power import generic_power

#import benchExpMCN


#### GESTION DU GRAPHE  - INITIALISATION ####

#file_AF = ".\Old_test\AF5_3.txt"
#file_AF = "./DATA/DAG_50_80/DAG_50_80_10.txt"
#file_AF = "./DATA/DAG_50_90/DAG_50_90_8.txt"
#file_AF = "./DATA/DAG_50_100/DAG_50_100_1.txt"

#file_AF = "./50_70_1500/DAG_50_70_760.txt"

#file_AF = "./DAG-test.txt"
#file_AF = "./small-test.txt"
#file_AF = "./delete-test.txt"
#file_AF = "./dep-test.txt"
#file_AF = "./50_75_1000/DAG_50_75_2.txt"

#file_AF = "./test-paper.txt"
#file_AF = "./diamond.txt"

#file_AF = "./DAG-Test-small.txt"      
#file_AF = "./diamond1.txt"     

file_AF = "./test1.txt"
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

def build_dico_paths(goal, ldico_att_in, dico_att):
    list_paths = get_paths_back(goal, ldico_att_in, dico_att, paths=None, current_path=None)
    #print(f"list paths: {list_paths}")
    print(f"nb paths: {len(list_paths)}")
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
            #if isinstance(dico_deg[att], sym.Basic)
            if type(dico_deg[att]) == sage.symbolic.expression.Expression or att in dep:
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


def nbTermes(liste_lvl, dico_att_in, dico_deg, dep):
    dico_term = {}
    dico_term_symb = {}
    i = len(liste_lvl)-1
    list_term = []
    som_terms = 0

    print(f"dep {dico_deg} and len {len(dico_deg)} ")

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
        som_terms = terms
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
                #t = t * 2 
            som_terms = som_terms + t

            list_term.append(t)
        
        print(f"list terms = {list_term}")
        print(f"max terms = {max(list_term)}")
        print(f"sum terms = {som_terms}\n")
        print(f"dico terms symb = {dico_term_symb}\n")

    if len(list_term) == 0:
        list_term.append(0)

    return dico_term,list_term,max(list_term),som_terms

############## ############## ############## ############## ############## 

def fastMCN(goal,dico_Arg,dico_Att,d_lvl):
    ldico_att_in = list_dico_Att_in(dico_Arg, dico_Att)
    dico_att_in = dlist_Att_to_Arg(ldico_att_in)
    dico_lvl = level_Arg(goal, 1, dico_att_in, d_lvl)
    liste_lvl = order_level(goal,dico_lvl)
    dico_deg = {}
    i = len(liste_lvl)-1
    #dico_paths = build_dico_paths(goal, ldico_att_in, dico_Att)

    #####
    dico_paths2, back_att = build_dico_paths2(goal, ldico_att_in, dico_Att)

    dico_att_out = list_dico_Att_out(dico_Arg, dico_Att)
    d_att_out =  dlist_Att_out_to_Arg(dico_att_out)
    
    #merge_dep = dependant_arg(goal, dico_att_in,d_att_out, dico_paths)

    merge_dep = dependant_arg2(d_att_out, back_att)

    print(f"dep = {merge_dep}")
    #print(f"dep2 = {merge_dep2}")
    
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

                string = str(dico_deg[goal])
                c1 = 0
                for elem2 in string:
                    if elem2 == "+" or elem2 == "-":
                        c1 += 1
                print(f"nb term avant à {arg} = {c1} ")

                #print(f"avant expr = {dico_deg[goal] }")

                #Expand
                start = time.time()
                dico_deg[goal] = dico_deg[goal].expand()
                end = time.time()
                elapsed = end - start

                #print(f"après expr = {dico_deg[goal] }")

                print(f"{arg} exp : {elapsed}")

                string = str(dico_deg[goal])
                c2 = 0
                for elem2 in string:
                    if elem2 == "+" or elem2 == "-":
                        c2 += 1
                print(f"nb term après à {arg} = {c2} ")
                print(f"diff = {c2-c1} ")

                #Delete
                for argu,lvl in liste_lvl:

                    for i in range(2,dico_deg[goal].degree(var(str(argu)))+1):
                        dico_deg[goal] = dico_deg[goal].subs(generic_power(var(str(argu)), i)==var(str(argu)))

                #Replace 
                dico_deg[goal] = dico_deg[goal].subs(var(str(arg))==dico_deg[arg]) 

    dico_term, list_term, max_list_term, som_terms = nbTermes(liste, dico_att_in, dico_deg, merge_dep)

    print(f"dico term {dico_term}")
    print(f"list term {list_term}")
    print(f"max_list_term {max_list_term}")
    print(f"som_terms {som_terms}")

    #conjunction_with_larger(goal, dico_att_in, dico_deg, merge_dep)

    return dico_deg[goal]

############################################################


start = time.time()
p = fastMCN("a",dico_Arg,dico_Att,dico_lvl)
end = time.time()
elapsed = end - start

print(f"Probability of a (ALL MCN algo) = {p}")
print(f"Time (ALL MCN algo) = {elapsed}\n")


"""

f = var('f')
e = var('e')
d = var('d')
c = var('c')

a = (1 - 0.2*c)*(1 - 0.1*(1 - 0.4*d)*(1 - 0.3*c))

a = a.expand()

print(f"first expand = {a}")

a2 = 0.9 - 0.156*(1 - 0.5*(0.2))*(1 - 0.6*d) + 0.04*d - 0.0176*(1 - 0.5*(0.2))*(1 - 0.6*d)*d

a2 = a2.expand()

print(f"second expand = {a2}")
"""