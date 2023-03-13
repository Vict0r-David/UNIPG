import time
import random
import numpy  
from itertools import combinations
import grounded as g


#### GESTION DU GRAPHE  - INITIALISATION ####

#file_AF = "AF_test2.txt"
#file_AF = ".\DAG_45\DAG_45_0.06_1.txt"
file_AF = ".\SCN_800\SCN_800_1.txt"
#file_AF = ".\Old_test\AF5_33.txt"
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


############# Approx #################

def one_Trial(arg,dico_Att,dico_Arg):
    count = 0
    random_att = {}
    for id, info in dico_Att.items():
        r = random.uniform(0,1)
        if -(info[2]) >= r:
            random_att[id] = info
    #if str(arg) in ground(random_att,dico_Arg):
    if str(arg) in g.grounded(dico_Arg, random_att):
        count += 1
    return count


def approx(arg,time_max,dico_Att,dico_Arg):
    count = 0
    n = 0
    sum_time = 0
    while sum_time < time_max:
        start = time.time()
        count += one_Trial(arg,dico_Att,dico_Arg)
        end = time.time()
        sum_time += end - start
        n += 1
    print(n)
    return float(count/n)


print(approx("a1",1,dico_Att,dico_Arg))


"""
############### Grounded ###########################

def ground(world,dico_Arg):
    d_arg = {}
    for arg in dico_Arg:
        d_arg[arg] = []
    
    In = []
    if len(world) > 0:
        for id,att in world.items():
            d_arg[att[1]].append(att[0])
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


def list_dico_Att_in(dico_Arg, dico_Att):
    ldico_att_in = {}
    for id_arg in dico_Arg:
        ldico_att_in[id_arg] = []
    for id,att in dico_Att.items():
        ldico_att_in[att[1]] += [id]
    return ldico_att_in


def dlist_Att_in_to_Arg(dico_list_att):
    dico = {}
    for arg,l in dico_list_att.items():
        dico[arg] = []
        for att in l:
            l2 = att.split("-")
            dico[arg].append(l2[0])
    return dico


def list_dico_Att_out(dico_Arg, dico_Att):
    ldico_att_out = {}
    for id_arg in dico_Arg:
        ldico_att_out[id_arg] = []
    for id,att in dico_Att.items():
        ldico_att_out[att[0]] += [id]
    return ldico_att_out

def dlist_Att_out_to_Arg(dico_list_att):
    dico = {}
    for arg,l in dico_list_att.items():
        dico[arg] = []
        for att in l:
            l2 = att.split("->")
            dico[arg].append(l2[1])
    return dico

def grounded(dico_Arg,dico_Att):  
    d_in = list_dico_Att_in(dico_Arg, dico_Att)
    att_in = dlist_Att_in_to_Arg(d_in)
    d_out = list_dico_Att_out(dico_Arg, dico_Att)
    att_out = dlist_Att_out_to_Arg(d_out)
    to_be_in = []
    label = {}
    und_pre = {}
    for x in dico_Arg:
        label[x] = "und"
        und_pre[x] = len(att_in[x]) 
        if und_pre[x] == 0:
            to_be_in.append(x)
 
    while len(to_be_in) > 0:
        x = to_be_in[0]
        to_be_in.remove(x)
        label[x] = "in"
        for y in att_out[x]:
            if label[y] != "out":
                label[y] = "out"
                for z in att_out[y]:
                    if label[z] == "und":
                        und_pre[z] = und_pre[z] - 1
                        if und_pre[z] == 0:
                            to_be_in.append(z)
    output = [] 
    for arg in label:
        if label[arg] == "in":
            output.append(arg)
    return output
"""
