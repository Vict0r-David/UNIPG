# Victor DAVID
# Created 03/10/2022

import math
import numpy  
import sympy as sym
from sympy import Lambda
from sympy import Pow
from sympy import arg
import time


#### GESTION DU GRAPHE  - INITIALISATION ####

#file_AF = ".\Old_test\AF5_3.txt"
#file_AF = "./DATA/DAG_50_80/DAG_50_80_10.txt"
#file_AF = "./DATA/DAG_50_90/DAG_50_90_10.txt"
#file_AF = "./DATA/DAG_50_100/DAG_50_100_1.txt"

#file_AF = "./DAG-test.txt"
file_AF = "./small-test.txt"
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

# def dependant_arg(goal,dico_att_in,dico_att_out,dico_paths):
#     dico_dep = {}
#     all_nodes = set()
#     for k in range(0,len(dico_att_in[goal])):
#         all_nodes.update(dico_paths[dico_att_in[goal][k]])
#     all_nodes.add(goal)
#     for att in dico_att_in[goal]:
#         dico_dep[att] = set()
#         branch = dico_paths[att]
#         i = 0
#         all_other_branchs = set()
#         while i < len(dico_att_in[goal]):
#             if dico_att_in[goal][i] != att:
#                 all_other_branchs.update(dico_paths[dico_att_in[goal][i]])
#             i += 1
#         inter = branch.intersection(all_other_branchs)
#         for arg in inter:
#             if len(dico_att_out[arg]) > 1:
#                 cpt = 0
#                 j = 0
#                 while cpt < 2 and j < len(dico_att_out[arg]):
#                     if dico_att_out[arg][j] in all_nodes:
#                         cpt += 1
#                     j += 1
#                 if cpt == 2:
#                         dico_dep[att].add(arg)
#     return dico_dep 

def dependant_arg(goal,dico_att_in,dico_att_out,dico_paths):
    #print(f"ICI : {dico_paths}")
    back_att = set()
    for arg,path in dico_paths.items():
        back_att.update(path)
    dep = set()
    #print(f"back att = {back_att}")
    for a in back_att:
        cpt = 0
        for b,path in dico_paths.items():
            if a in path:
                cpt += 1
            if cpt > 1:
                cpt2 = 0
                for c in dico_att_out[a]:
                    #print(f"for {a} we test {c}")
                    if c in back_att:
                        cpt2 += 1
                    if cpt2 > 1:
                        #print(f"add {a}")
                        dep.add(a)
                        break
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

def nbTermes(liste_lvl, dico_att_in, dico_deg, dep):
    dico_term = {}
    i = len(liste_lvl)-1

    while i >= 0: 
        l_att = dico_att_in[str(liste_lvl[i][0])]
        dico_term[str(liste_lvl[i][0])] = 1
        if isinstance(dico_deg[str(liste_lvl[i][0])], sym.Basic):
            for att in l_att:
                if isinstance(dico_deg[att], sym.Basic) or (att in dep):
                    dico_term[str(liste_lvl[i][0])] = dico_term[str(liste_lvl[i][0])] * (dico_term[att]+1)
        i-=1
    return dico_term 


def nbTermes2(liste_lvl, dico_att_in, dico_deg, dep, taille):
    dico_term = {}
    i = len(liste_lvl)-1

    while i >= 0: 
        l_att = dico_att_in[str(liste_lvl[i][0])]
        dico_term[str(liste_lvl[i][0])] = 1
        if isinstance(dico_deg[str(liste_lvl[i][0])], sym.Basic):
            for att in l_att:
                if att in dep:
                    dico_term[str(liste_lvl[i][0])] = dico_term[str(liste_lvl[i][0])] * 2
                elif isinstance(dico_deg[att], sym.Basic):
                    dico_term[str(liste_lvl[i][0])] = dico_term[str(liste_lvl[i][0])] * (dico_term[att]+1)
        i-=1
    

    if len(dep) > 0:
        j = 0
        terms =  dico_term[str(liste_lvl[0][0])]
        # keep only the dependant argument and thanks to list_lvl they are ordered increasingly
        while j < len(liste_lvl):
            if liste_lvl[j][0] not in dep:
                del liste_lvl[j]
                j -=1
            j +=1
        #print(liste_lvl)
        minus = 0
        for arg,lvl in liste_lvl:   
            terms = terms + taille*2**(len(liste_lvl)-minus) + 2**(len(liste_lvl)-minus-1) * (dico_term[arg]-1)

        print(f"nb terms = {terms}")

    return dico_term 


def symb_arg(goal,dico_paths,dep,d_att_out):
    dico_symb = {}
    for arg in dep:
        dico_symb[arg] = set()
    for a in dico_paths[goal]:
        for arg in dep:
            if arg in dico_paths[a] and a != goal and a not in dep :
                dico_symb[arg].add(a)
    return dico_symb

def size_symb(arg,symb_arg):
    cpt = 0
    for dep,liste in symb_arg.items():
        if arg in liste:
            cpt += 1
    return cpt

def approx_comlexity(goal,set_symb,dico_symb,dep):
    val = 1
    for symb in set_symb:
        #val = val * (size_symb(symb,dico_symb) + 1)
        val = val * (2)
    for d in dep:
        val = val * 2
    val = val - 1
    if 2**(len(dep)-1) > 0: 
        val = val * 2**(len(dep)-1)
    val += 2**(len(dep))
    return val

def dep_att(dep,d_att_out):
    set_dep = set()
    for arg in dep:
        set_dep.add((arg,len(d_att_out[arg])))
    return set_dep

############## ############## ############## ############## ############## 

def fastMCN(goal,dico_Arg,dico_Att,d_lvl):
    ldico_att_in = list_dico_Att_in(dico_Arg, dico_Att)
    dico_att_in = dlist_Att_to_Arg(ldico_att_in)
    dico_lvl = level_Arg(goal, 1, dico_att_in, d_lvl)
    liste_lvl = order_level(goal,dico_lvl)
    dico_deg = {}
    i = len(liste_lvl)-1
    dico_paths = build_dico_paths(goal, ldico_att_in, dico_Att)
    dico_att_out = list_dico_Att_out(dico_Arg, dico_Att)
    d_att_out =  dlist_Att_out_to_Arg(dico_att_out)
    # dico_dep = dependant_arg(goal, dico_att_in,d_att_out, dico_paths)
    # merge_dep = set()
    # for id,setdep in dico_dep.items():
    #     merge_dep.update(setdep)
    
    merge_dep = dependant_arg(goal, dico_att_in,d_att_out, dico_paths)

    #dico_symb = symb_arg(goal,dico_paths,merge_dep,d_att_out)
    #print(f"dico symb: {dico_symb}")
    #set_dep = dep_att(merge_dep,d_att_out)
    #print(f"set dep: {set_dep}")
    
    symbo = False
    while i >= 0: 
        l_att = ldico_att_in[str(liste_lvl[i][0])]
        deg = 1
        for att in l_att:
                if dico_Att[att][0] in merge_dep :
                    deg = deg * (1-(sym.Symbol(dico_Att[att][0])* -dico_Att[att][2]))
                    symbo = True
                else:
                    deg = deg * (1-(dico_deg[dico_Att[att][0]]*- dico_Att[att][2])) 
        dico_deg[str(liste_lvl[i][0])] = deg
        i-=1
    if symbo:
        dico_deg[goal] = sym.expand(dico_deg[goal])
        dico_deg[goal] = dico_deg[goal].replace(Pow, lambda a,b: Pow(a,1))

    
    # set_symb = set()
    # for arg,symbo in dico_symb.items():
    #     set_symb.update(symbo)

    #val = approx_comlexity(goal,set_symb,dico_symb,set_dep)
    #print(f"approx complexity = {val}")

    #print(f"nb symb = {len(set_symb)}")
    
    
    # moy_symb = 0
    # max_symb = 0
    # it = 1
    # a_in = 1
    # for arg,liste in dico_symb.items():
    #     l = len(liste)
    #     it += 1
    #     if l > max_symb:
    #         max_symb = l
    #     moy_symb += l
    #     if len(dico_att_in[arg]) > 1:
    #         a_in = a_in * len(dico_att_in[arg])
    # moy_symb = moy_symb/it
    # print(f"moy nb arg symb = {moy_symb}")
    # print(f"max nb symb = {max_symb}")
    # print(f"prod att in symb = {a_in}")

    print(f"nb dep = {len(merge_dep)}")
    print(f"dep = {merge_dep}")
    # moy_dep = 0
    # max_dep = 0
    # it = 1
    # a_out = 1
    # for tup in set_dep:
    #     nbAtt = tup[1]
    #     it += 1
    #     if nbAtt > max_dep:
    #         max_dep = nbAtt
    #     moy_dep += nbAtt
    #     a_out = a_out * nbAtt
    # moy_dep = moy_dep/it
    # print(f"moy nb att dep = {moy_dep}")
    # print(f"max nb att dep = {max_dep}")
    # print(f"prod att out dep = {a_out}")

    dico_conj = conjunction_with_larger(goal, dico_att_in, dico_deg, merge_dep)
    print(f"dico conjonction = {dico_conj}\n")

    #dico_terme = nbTermes(liste_lvl, dico_att_in, dico_deg, merge_dep)
    #print(f"dico termes = {dico_terme}")

    liste = list(liste_lvl)

    if len(merge_dep) > 0:
            j = 0
            # keep only the dependant argument and thanks to list_lvl they are ordered increasingly
            while j < len(liste_lvl):
                if liste_lvl[j][0] not in merge_dep:
                    del liste_lvl[j]
                    j -=1
                j +=1
            #print(liste_lvl)
            for arg,lvl in liste_lvl:   
                dico_deg[goal] = dico_deg[goal].subs(arg,dico_deg[arg]) 
                dico_deg[goal] = sym.expand(dico_deg[goal])
                dico_deg[goal] = dico_deg[goal].replace(Pow, lambda a,b: Pow(a,1))


    taille = len(str(dico_deg[goal]))
    dico_terme2 = nbTermes2(liste, dico_att_in, dico_deg, merge_dep, 1)
    print(f"dico termes  = {dico_terme2}")

    return dico_deg[goal]

############################################################


start = time.time()
#a45
p = fastMCN("a",dico_Arg,dico_Att,dico_lvl)
end = time.time()
elapsed = end - start

print(f"Probability of a (ALL MCN algo) = {p}")
print(f"Time (ALL MCN algo) = {elapsed}\n")



"""
x = sym.Symbol('x')
y = sym.Symbol('y')
a = 4 + x**2 + y
b = 4 + 5**2

print(isinstance(a, sym.Basic))
print(isinstance(b, sym.Basic))



da = sym.Symbol("da")
ab = sym.Symbol("ab")
d = sym.Symbol("d")
exp = (1 - (1 - d * da) * ab) * (1 - d * da)
exp = exp.expand()
print(exp)
"""