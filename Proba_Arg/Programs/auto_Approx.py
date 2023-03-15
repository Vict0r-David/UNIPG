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
    #print(n)
    return float(count/n),n


#print(approx("a1",1,dico_Att,dico_Arg))



def init_dico_lvl(dico_Arg,dico_lvl):
    for arg in dico_Arg:
        dico_lvl[arg] = 0
    return dico_lvl


#### GESTION DU GRAPHE  - INITIALISATION ####

output = ".\DAG_20+\output_Approx_DAG+_20.txt"
#output = ".\SCN_200\output_UF_SCN_200_1.txt"
f_out = open(output,"w")
avg_total_time = 0
avg_total_nodes = 0
avg_total_edges = 0
total_max = 0
total_it = 0

for i in range(1,21):
    print(i)
    file_AF = ".\DAG_20+\DAG_20_0.2_"+str(i)+".txt"
    #file_AF = ".\SCN_200\SCN_200_"+str(i)+".txt"
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
    avg_it = 0
    for a in dico_Arg:
        print(a)
        start = time.time()
        #proba = Fast(str(a), dico_Arg, dico_Att, dico_lvl, dico_pw)
        proba,it = approx(str(a),0.78,dico_Att,dico_Arg)
        #proba = AllfastMCN(str(a),dico_Arg,dico_Att,dico_lvl)
        end = time.time()
        avg_it += it
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
    f_out.write("Average_Iterations = ")
    avg_it = avg_it/len(liste_output)
    f_out.write(str(avg_it))
    nb_nodes = len(dico_Arg)
    nb_edges = len(dico_Att)
    f_out.write("\nnb_nodes = ")
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
    total_it += avg_it
    AF.close()

avg_total_time = avg_total_time/20
avg_total_nodes = avg_total_nodes/20
avg_total_edges = avg_total_edges/20
total_it = total_it/20
f_out.write("Average Total Time = ")
f_out.write(str(avg_total_time))
f_out.write(" Average Total Iterations = ")
f_out.write(str(total_it))
f_out.write("\nAverage Total Nodes = ")
f_out.write(str(avg_total_nodes))
f_out.write("\nAverage Total Edges = ")
f_out.write(str(avg_total_edges))
f_out.write("\nMax Total Time for one argument = ")
f_out.write(str(total_max))
f_out.close()