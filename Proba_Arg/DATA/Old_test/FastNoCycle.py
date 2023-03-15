# Victor DAVID
# Created 03/10/2022

import math
import numpy  
import sympy as sym

def grad(dico_Arg,dico_Att):
    new_arg = {}
    for id_a,w_ar in dico_Arg.items():
        deg = 1
        for id_att,att in dico_Att.items():
            if id_a == att[1]:
               deg = deg * (1-(numpy.float64(dico_Arg[att[0]])*numpy.float64(-att[2])))
        new_arg[id_a] = deg
    return new_arg

def test_end(dico_Arg, dico_Arg_pred, threshold):
    for arg,w_ar in dico_Arg.items():
        if math.sqrt((dico_Arg_pred[arg]- w_ar)**2 + (w_ar - dico_Arg_pred[arg])**2) > threshold:
            return False
    return True

def semantics(dico_Arg, dico_Att, threshold):
    dico_Arg_pred = dico_Arg
    dico_Arg = grad(dico_Arg, dico_Att)
    
    while not test_end(dico_Arg, dico_Arg_pred, threshold):
        print(dico_Arg)
        dico_Arg_pred = dico_Arg
        dico_Arg = grad(dico_Arg, dico_Att)
    return dico_Arg


file_AF = "AF_0.txt"
AF = open(file_AF,"r")
dico_Arg = {}
dico_Att = {}
for line in AF:
    if line[0:3] == "arg":
        l1 = line.partition("(")
        l2 = l1[2].partition(")")
        dico_Arg[l2[0]] = 1
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

threshold = 0.000001
print(semantics(dico_Arg, dico_Att, threshold))



################################################################################################################
############################################### Usefull Function ###############################################

def list_dico_Att_out(dico_Arg, dico_Att):
    ldico_att_out = {}
    for id_arg in dico_Arg:
        ldico_att_out[id_arg] = []
    for id,att in dico_Att.items():
        ldico_att_out[att[0]] += [id]
    return ldico_att_out

def list_dico_Att_in(dico_Arg, dico_Att):
    ldico_att_in = {}
    for id_arg in dico_Arg:
        ldico_att_in[id_arg] = []
    for id,att in dico_Att.items():
        ldico_att_in[att[1]] += [id]
    return ldico_att_in


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


def psum(l_value):
    if len(l_value) == 0:
        return 0
    elif len(l_value) == 1:
        return l_value[0]
    else:
        x = (l_value[0]+ l_value[1]) - (l_value[0] * l_value[1])
        for i in range(2,len(l_value)):
            x = (x+ l_value[i]) - (x * l_value[i])
        return x

def prod(list_value):
    p = 1
    for v in list_value:
        p = p*v
    return p


################################################################################################################
################################################################################################################

############################################### Controversial Case ##############################################


# Symbolic Computation 

x = sym.Symbol('x')
y = sym.Symbol('y')
res = sym.expand((x+y)*x * y**10 - x**2)
stres = str(res)
print(stres)

res2 = sym.expand(1 - (1 - x*0.4) * 0.8 )
stres2 = str(res2)
print(stres2)
#l_res = stres.split("+")

def reducePower(string):
    entier = ['0','1','2','3','4','5','6','7','8','9']
    while "**" in string:        
        i = string.index("**")
        string = string[:i] + string[i+2:]
        while i < len(string) and string[i] in entier:
            if i+1<len(string):
                string = string[:i] + string[i+1:]
            else:
                string = string[:i]
    return string

print(reducePower(stres))

############################################################

# Preprocessing

#For any level, i.e. any node
def N_PathNumber(node, dico_Arg, dico_Att):
    l_att_in = list_dico_Att_in(dico_Arg, dico_Att)
    path_back = get_paths_back(node, l_att_in, dico_Att)
    list_dico_forall_path = []
    set_arg = set()

    for path in path_back:
        dico_temp = {}
        for arg in path:
            set_arg.add(arg)
            dico_temp[arg] = []
            for id in dico_temp:
                if arg != id:
                    dico_temp[id].append([arg,len(dico_temp[id])+1])
        list_dico_forall_path.append(dico_temp)
    
    return (list_dico_forall_path, set_arg, path_back)


def Find_un_controversial(node, l_path, set_arg, path_back, dico_Att):    
    attDef_dico = {}
    for arg in set_arg:
        attDef_dico[arg] = [0,0]
    for path in l_path:
        for pair in path:
            if pair[1]%2 == 0:
                attDef_dico[pair[0]][0] += 1
            else:
                attDef_dico[pair[0]][1] += 1
    Reas = set()
    idUnreas = set()
    Unreas = {}
    for arg in set_arg:
        if attDef_dico[arg][0] > 0 and attDef_dico[arg][1] > 0:
            idUnreas.add(arg)
            Unreas[arg] = set()
        elif attDef_dico[arg][0] == 0 and attDef_dico[arg][1] > 0:
            id = str(arg) + "->" + str(node)
            if id in dico_Att:
                Reas.add(arg)
    for path in path_back:
        set_temp = set()
        for arg in path:
            set_temp.add(arg)
            if arg in idUnreas:
                Unreas[arg].update(set_temp)
    return Reas,Unreas

def N_OrderUnReas(node, dico_Arg, dico_Att):
    list_dico_path, set_arg, path_back = N_PathNumber(node, dico_Arg, dico_Att)
    dico_output = {}
    dico_path = {}
    max_dico = {}
    dico_order = {}
    l_arg = list(set_arg)
    for i in range(0,len(set_arg)):
        max_dico[l_arg[i]] = 0
        dico_path[l_arg[i]] = []
        dico_order[i] = []

    for d_path in list_dico_path:
        for arg,liste in d_path.items():
            if len(liste) > 0:
                dico_path[arg].append(liste)

    for arg,l_path in dico_path.items():
        for path in l_path:
            for pair in path:
                if pair[1] > max_dico[pair[0]]:
                    max_dico[pair[0]] = pair[1]     

    for arg in set_arg:
        reas,unreas = Find_un_controversial(arg,dico_path[arg],set_arg,path_back,dico_Att)
        dico_output[arg] = [reas,unreas]
        dico_order[max_dico[arg]].append(arg)       

    return dico_order,dico_output









"""
# For node, i.e. last level
def pathNumber(node, dico_Arg, dico_Att):
    l_att_in = list_dico_Att_in(dico_Arg, dico_Att)
    path_back = get_paths_back(node, l_att_in, dico_Att)
    l_out = []
    set_arg = set()
    set_arg.add(node)
    for path in path_back:
        l_temp = []
        for i in range(1,len(path)):
            set_arg.add(path[i])
            l_temp.append([path[i],i])
        l_out.append(l_temp)
    return (l_out, set_arg, path_back)

def OrderUnReas(node, dico_Arg, dico_Att):
    l_path, set_arg, path_back = pathNumber(node, dico_Arg, dico_Att)
    
    max_dico = {}
    attDef_dico = {}
    for arg in set_arg:
        max_dico[arg] = 0
        attDef_dico[arg] = [0,0]
    for path in l_path:
        for pair in path:
            if pair[1] > max_dico[pair[0]]:
                max_dico[pair[0]] = pair[1] 
            if pair[1]%2 == 0:
                attDef_dico[pair[0]][0] += 1
            else:
                attDef_dico[pair[0]][1] += 1

    dico_order = {}
    for i in range(0,len(set_arg)):
        dico_order[i] = []
    Reas = set()
    idUnreas = set()
    UnReas = {}

    for arg in set_arg:
        dico_order[max_dico[arg]].append(arg)
        if attDef_dico[arg][0] > 0 and attDef_dico[arg][1] > 0:
            idUnreas.add(arg)
            UnReas[arg] = set()
        elif attDef_dico[arg][0] == 0 and attDef_dico[arg][1] > 0:
            id = str(arg) + "->" + str(node)
            if id in dico_Att:
                Reas.add(arg)
    
    for path in path_back:
        set_temp = set()
        for arg in path:
            set_temp.add(arg)
            if arg in idUnreas:
                UnReas[arg].update(set_temp)

    return dico_order,Reas,UnReas
"""

############################################################
"""
l_att_in = list_dico_Att_in(dico_Arg, dico_Att)
path_back = get_paths_back("a1", l_att_in, dico_Att)

for path in path_back:
    print(path)

dico_paths,set_arg, path_back = pathNumber("a1", dico_Arg, dico_Att)

print(set_arg)

for dico in dico_paths:
    print(dico)

dicoOrder,Reas,UnReas = OrderUnReas("a1", dico_Arg, dico_Att)
print(dicoOrder)
print(Reas)
print(UnReas)
"""

d_order, output = N_OrderUnReas("a1", dico_Arg, dico_Att)
print(f'Order: {d_order} \n[arg: Reas(arg), UnReas(arg)]')

for id,item in output.items():
    print(f'{id}: {item}')


