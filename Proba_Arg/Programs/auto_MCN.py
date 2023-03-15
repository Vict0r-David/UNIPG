# Victor DAVID
# Created 03/10/2022

import time
import numpy
#import Fast_MCN as mcn
#import Fast_SCN as scn
#import Approx as app
#import constellation as con

from sympy import Pow
import sympy as sym


###########################################################################################################


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
    dico_paths = {}
    for path in list_paths:
        for i in range(0,len(path)):
            if str(path[i]) not in dico_paths:
                dico_paths[path[i]] = set(path[i:])
            else:
                dico_paths[str(path[i])].update(set(path[i:]))
    return dico_paths

def dependant_arg(goal,dico_att_in,dico_att_out,dico_paths):
    dico_dep = {}
    all_nodes = set()
    for k in range(0,len(dico_att_in[goal])):
        all_nodes.update(dico_paths[dico_att_in[goal][k]])
    all_nodes.add(goal)
    for att in dico_att_in[goal]:
        dico_dep[att] = set()
        branch = dico_paths[att]
        i = 0
        all_other_branchs = set()
        while i < len(dico_att_in[goal]):
            if dico_att_in[goal][i] != att:
                all_other_branchs.update(dico_paths[dico_att_in[goal][i]])
            i += 1
        inter = branch.intersection(all_other_branchs)
        for arg in inter:
            if len(dico_att_out[arg]) > 1:
                cpt = 0
                j = 0
                while cpt < 2 and j < len(dico_att_out[arg]):
                    if dico_att_out[arg][j] in all_nodes:
                        cpt += 1
                    j += 1
                if cpt == 2:
                        dico_dep[att].add(arg)
    return dico_dep 

def fastMCN(goal,dico_Arg,dico_Att,d_lvl,dico_deg_computed):
    ldico_att_in = list_dico_Att_in(dico_Arg, dico_Att)
    dico_att_in = dlist_Att_to_Arg(ldico_att_in)
    dico_lvl = level_Arg(goal, 1, dico_att_in, d_lvl)
    liste_lvl = order_level(goal,dico_lvl)
    dico_deg = {}
    i = len(liste_lvl)-1
    dico_paths = build_dico_paths(goal, ldico_att_in, dico_Att)
    dico_att_out = list_dico_Att_out(dico_Arg, dico_Att)
    d_att_out =  dlist_Att_out_to_Arg(dico_att_out)
    dico_dep = dependant_arg(goal, dico_att_in,d_att_out, dico_paths)
    
    merge_dep = set()
    for id,setdep in dico_dep.items():
        merge_dep.update(setdep)
           
    symbo = False
    while i >= 0: 
        l_att = ldico_att_in[str(liste_lvl[i][0])]
        deg = 1
        for att in l_att:
                if dico_Att[att][0] in merge_dep :
                    deg = deg * (1-(sym.Symbol(dico_Att[att][0])* -dico_Att[att][2]))
                    symbo = True
                # if  the attacker have some variables, they are dependant and so we need to keep it     
                elif isinstance(dico_deg[dico_Att[att][0]], sym.Basic) and type(dico_deg[dico_Att[att][0]]) != sym.core.numbers.Float:
                    deg = deg * (1-(dico_deg[dico_Att[att][0]]*- dico_Att[att][2]))
                # otherwise the attacker do not have any variables, therefore we can replace the attacker by its final probability
                else:
                    deg = deg * (1-(dico_deg_computed[dico_Att[att][0]]*- dico_Att[att][2])) 
        dico_deg[str(liste_lvl[i][0])] = deg
        i-=1
    if symbo:
        dico_deg[goal] = sym.expand(dico_deg[goal])
        dico_deg[goal] = dico_deg[goal].replace(Pow, lambda a,b: Pow(a,1))
        

    return dico_deg[goal], merge_dep, liste_lvl

def AllfastMCN(goal,dico_Arg,dico_Att,dico_lvl):
    ldico_att_in = list_dico_Att_in(dico_Arg, dico_Att)
    dico_att_in = dlist_Att_to_Arg(ldico_att_in)
    dico_lvl_init = level_Arg(goal, 1, dico_att_in, dico_lvl)
    liste_lvl = order_level(goal,dico_lvl_init)
    dico_deg_computed = {}
    for arg in dico_Arg:
        dico_deg_computed[arg] = 1
    i = len(liste_lvl)-1

    while i >= 0: 
        deg, set_dep, liste_lvl_j = fastMCN(liste_lvl[i][0],dico_Arg,dico_Att,dico_lvl,dico_deg_computed)
        dico_deg_computed[liste_lvl[i][0]] = deg
        if len(set_dep) > 0:
            j = 0
            # keep only the dependant argument and thanks to list_lvl they are ordered increasingly
            while j < len(liste_lvl_j):
                if liste_lvl_j[j][0] not in set_dep:
                    del liste_lvl_j[j]
                    j -=1
                j +=1
            
            for arg,lvl in liste_lvl_j:   
                # because if during the compution deg = 0, it is useless to expand more and it is not possible to do deg.subs on int           
                if type(deg) != int :           
                    new_arg = 1
                    if len(ldico_att_in[arg]) > 0:
                        for att in ldico_att_in[arg]:
                            new_arg = new_arg * (1- sym.Symbol(dico_Att[att][0]) * float(-dico_Att[att][2]))
                        deg = deg.subs(arg,new_arg) 
                        sum_term = 0
                        for term in deg.args :
                            for var in term.free_symbols:
                                if str(var) not in set_dep or dico_deg_computed[str(var)] == 1:
                                    term = term.subs(var,dico_deg_computed[str(var)])
                            term = sym.expand(term)
                            term = term.replace(Pow, lambda a,b: Pow(a,1))
                            sum_term = sum_term + term
                        deg = sum_term
                    else :
                        deg = deg.replace(sym.Symbol(arg),1)

            dico_deg_computed[liste_lvl[i][0]] = deg
        i -= 1
    return dico_deg_computed[liste_lvl[0][0]]



###########################################################################################################
############################################# Experimentation #############################################

def init_dico_lvl(dico_Arg,dico_lvl):
    for arg in dico_Arg:
        dico_lvl[arg] = 0
    return dico_lvl


#### GESTION DU GRAPHE  - INITIALISATION ####

#output = ".\DAG_50_100\output_MCN_DAG_50_100.txt"
output = ".\SCN_1000\output_M_SCN_1000.txt"
f_out = open(output,"w")
avg_total_time = 0
avg_total_nodes = 0
avg_total_edges = 0
total_max = 0

for i in range(1,2):
    print("")
    print(i)
    #file_AF = ".\DAG_50_100\DAG_50_100_"+str(i)+".txt"
    file_AF = ".\SCN_1000\SCN_1000_"+str(i)+".txt"
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

    max_time = 0
    min_time = 999999999999999999
    liste_output = []
    for a in dico_Arg:
        print(a)
        start = time.time()
        proba = AllfastMCN(str(a),dico_Arg,dico_Att,dico_lvl)
        end = time.time()
        elapsed = end - start
        liste_output.append([proba,round(elapsed,4)])
        dico_lvl = init_dico_lvl(dico_Arg,dico_lvl)

        if elapsed > max_time:
            max_time = elapsed
        if elapsed < min_time:
            min_time = elapsed
        
        #print(proba,round(elapsed,4))

    #print(f'Probability of a = {proba}')
    #print(f'Time: {elapsed:.5}s\n')

        f_out.write("Prob arg(")
        f_out.write(str(a))
        f_out.write(") = ")
        f_out.write(str(proba))
        f_out.write(" and Time = ")
        f_out.write(str(round(elapsed,4)))
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

avg_total_time = avg_total_time/10
avg_total_nodes = avg_total_nodes/10
avg_total_edges = avg_total_edges/10
f_out.write("Average Total Time = ")
f_out.write(str(avg_total_time))
f_out.write("\nAverage Total Nodes = ")
f_out.write(str(avg_total_nodes))
f_out.write("\nAverage Total Edges = ")
f_out.write(str(avg_total_edges))
f_out.write("\nMax Total Time for one argument = ")
f_out.write(str(total_max))
f_out.close()