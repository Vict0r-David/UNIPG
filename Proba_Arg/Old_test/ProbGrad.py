# Victor DAVID
# Created 03/10/2022

import math
import numpy  
import time 

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
        #print(dico_Arg)
        dico_Arg_pred = dico_Arg
        dico_Arg = grad(dico_Arg, dico_Att)
    return dico_Arg


#file_AF = "AF5_4.txt"
file_AF = ".\SCN_1000\SCN_1000_1.txt"
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

threshold = 0.0001
start = time.time()
dico = semantics(dico_Arg, dico_Att, threshold)
end = time.time()
elapsed = end - start
#print(dico["a1"])


print(f"Probability of a1 = {dico['a1']}")
print(f"Time  = {elapsed}\n")
