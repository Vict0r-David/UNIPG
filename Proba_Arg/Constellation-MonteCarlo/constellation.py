
# Created 20/10/2022

import math 
import numpy
#from itertools import combinations
from itertools import chain, combinations, permutations
import time
import grounded as g

def powerset(iterable):
    list_combinations = list()
    for n in range(len(iterable) + 1):
        list_combinations += list(combinations(iterable, n))
    return list_combinations


def ground(world,dico_Arg):
    d_arg = {}
    for arg in dico_Arg:
        d_arg[arg] = []
    
    In = []

    if len(world) > 0:
        for att in world:
            d_arg[att[1][1]].append(att[1][0])
        #print(d_arg)
        Out = [] 
        Stop = False
        while len(In)+len(Out) < len(d_arg) and not Stop:
            nb = len(In)+len(Out)
            for arg,att in d_arg.items():
                if arg not in In and arg not in Out:
                    goIn = True
                    for attack in att:
                        if attack not in Out:
                            goIn = False
                            break
                    if goIn:
                        In.append(arg)
                    else:
                        goOut = False
                        for attack in att:
                            if attack in In:
                                goOut = True
                        if goOut:
                            Out.append(arg) 
            nb2 = len(In)+len(Out)
            if nb2-nb == 0:
                Stop = True
    else:
        for arg in dico_Arg:
            In.append(arg)
    return In

def world_to_attack(world):
    dico = {}
    for elem in world:
        dico[elem[0]] = elem[1]
    return dico

def build_Constellation(dico_Att,dico_Arg):
    set_att = set()
    set_attinfo = set()
    for att,info in dico_Att.items():
        attinfo = [att,tuple(info)]
        attinfo = tuple(attinfo)
        set_att.add(att)
        set_attinfo.add(attinfo)
    list_worlds = powerset(set_attinfo)

    constellation = []
    for world in list_worlds:
        proba = 1
        seen_att = set()
        for l_att in world:
            seen_att.add(l_att[0])
            proba = numpy.float64(proba) * numpy.float64(-l_att[1][2])
        for att in (set_att - seen_att):
            proba = numpy.float64(proba) * numpy.float64(1 + dico_Att[att][2])

        In = ground(world,dico_Arg)
        #d_att = world_to_attack(world)
        #In = g.grounded(dico_Arg,d_att)

        constellation.append([proba,In,world])
    return constellation

def proba_arg(dico_Att,dico_Arg):
    d_arg = {}
    for arg in dico_Arg:
        d_arg[arg] = 0
    constellation = build_Constellation(dico_Att,dico_Arg)
    for world in constellation:
        for arg in d_arg:
            if arg in world[1]:
                d_arg[arg] += world[0]
    return d_arg


#file_AF = "AF_0.txt"
#file_AF = ".\Old_test\AF5_3.txt"
#file_AF = ".\DAG_50_21\DAG_50_21_3.txt"
#file_AF = "./././verif.txt"
file_AF = "./test1.txt"
#file_AF = "./small-test.txt"
#file_AF = "./delete-test.txt"
#file_AF = "./dep-test.txt"

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

#threshold = 0.0001

solution = []

Const = build_Constellation(dico_Att,dico_Arg)
for w in Const:
    list_att = []
    for att in w[2]:
        list_att.append((att[0],att[1][2]))
    ###########if 'b' in w[1] and (('d12->b', -0.3) not in list_att and (('d2->b', -0.6) in list_att and ('a->b', -0.2) in list_att)):
    #print(f'{w[0]}, {w[1]}, {list_att}') 
    solution.append(round(w[0],12))
    

start = time.time()
dico_proba = proba_arg(dico_Att,dico_Arg)
end = time.time()
elapsed = end - start

print(dico_proba)
print(f'Time: {elapsed:.5}s\n')


#print(solution)

#print(len(solution))


def sum_weight(solution):
    sum = 0
    for v in solution:
        sum += v
    return sum

#print(sum_weight(solution))
