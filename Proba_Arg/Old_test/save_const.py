# Victor DAVID
# Created 20/10/2022

import math 
from itertools import combinations

def powerset(iterable):
    list_combinations = list()
    for n in range(len(iterable) + 1):
        list_combinations += list(combinations(iterable, n))
    return list_combinations



def build_Constellation(dico_Att):
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
            proba = proba * (-l_att[1][2])
        for att in (set_att - seen_att):
            proba = proba * (1 + dico_Att[att][2])
        constellation.append([round(proba,8),world])
    return constellation


def grounded(world,dico_Arg):
    d_arg = {}
    for arg in dico_Arg:
        d_arg[arg] = []
        
    for att in world[1]:
        d_arg[att[1][1]].append(att[1][0])
    In = []    
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
    return In



def proba_arg(dico_Att,dico_Arg):
    d_arg = {}
    for arg in dico_Arg:
        d_arg[arg] = 0
    constellation = build_Constellation(dico_Att)
    for world in constellation:
        In = grounded(world,dico_Arg)
        for arg in d_arg:
            if arg in In:
                d_arg[arg] += world[0]
    return d_arg


file_AF = "AF5_4.txt"
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
        w = float(l4[0][:-1])
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
Const = build_Constellation(dico_Att)
for c in Const:
    print(c)

dico_proba = proba_arg(dico_Att,dico_Arg)
print(dico_proba)