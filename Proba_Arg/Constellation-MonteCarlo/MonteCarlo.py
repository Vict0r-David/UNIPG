import time
import random
import numpy  
from itertools import combinations
import grounded as g


#### GESTION DU GRAPHE  - INITIALISATION ####

#file_AF = "AF_test2.txt"
#file_AF = ".\DAG_20+\DAG_20_0.2_12.txt"
#file_AF = ".\SCN_1000\SCN_1000_3.txt"
#file_AF = ".\Old_test\AF5_33.txt"
"""
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
"""

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
    #print(n)
    return float(count/n)


#print(approx("a1",0.071, dico_Att,dico_Arg))
