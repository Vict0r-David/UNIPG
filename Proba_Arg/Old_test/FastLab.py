# Victor DAVID
# Created 03/10/2022

import math
import numpy  

def grad(dico_Arg,dico_Att):
    new_arg = {}
    dico_in = list_dico_Att_in(dico_Arg, dico_Att)
    print(dico_in)
    for id_a in dico_Arg:
        deg = 1
        for att in dico_in[id_a]:
            deg = deg * (1- (dico_Arg[dico_Att[att][0]] * (- dico_Att[att][2]) )) 
            #new_arg[id_a] = deg
        #new_arg[id_a] = deg - error(id_a, dico_Arg, dico_Arg, dico_Att)
        err = error(id_a, dico_Arg, dico_Arg, dico_Att)
        new_arg[id_a] = deg #- err
        if id_a == 'b':
            print(f'{err} \n')
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

def error(node, dico_deg, dico_Arg, dico_att):
    err = 0
    l_att_in = list_dico_Att_in(dico_Arg, dico_att)
    l_att_out = list_dico_Att_out(dico_Arg, dico_att)
    for id_att in l_att_in[node]:
        node_attack = dico_att[id_att][0]
        att_base = all_att(node_attack, node, l_att_out, dico_att)
        base = 1
        for id_att_base in att_base:
            base = base * (-dico_att[id_att_base][2])

        l_att, l_def = path_att_def(node_attack, node, l_att_out, dico_att)
        l_prob_att = []
        l_probInv_att = []
        for id_att_node in l_att_in[dico_att[id_att][0]]:
            l_prob_att.append( dico_deg[dico_att[id_att_node][0]] * (-dico_att[id_att_node][2]) )
            l_probInv_att.append(1- dico_deg[dico_att[id_att_node][0]] * (-dico_att[id_att_node][2]) )
        
        prob_att = psum(l_prob_att) ** len(l_att)
        if len(l_def) == 0:
            prob_def = 0
        else:
            prob_def = prod(l_probInv_att) ** len(l_def)

        err += base * prob_att * prob_def

        if node == 'b':
            print(f'{base} * {prob_att} * {prob_def}')

    return err



file_AF = "AF5_2.txt"
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

threshold = 0.00000001
print(semantics(dico_Arg, dico_Att, threshold))


#ldico_att = list_dico_Att_out(dico_Arg, dico_Att)
#print(ldico_att)

"""
paths = get_paths('e', 'b', ldico_att, dico_Att)
print(paths)

gp = good_paths('e', 'b', ldico_att, dico_Att)
print(gp)

all_a = all_att('e', 'b', ldico_att, dico_Att)
print(all_a)

l_att,l_def = path_att_def('e', 'b', ldico_att, dico_Att)
print(l_att)
print(l_def)

print(psum([0.7,0.8,0.9]))


dico_deg = {}
dico_deg['c'] = 1
dico_deg['d'] = 0.1
dico_deg['a'] = 0.93

#print(error('b', dico_deg, dico_Att))
"""