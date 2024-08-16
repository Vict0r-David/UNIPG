import math 
import pickle
import numpy

import MonteCarlo as MC

import matplotlib.pyplot as plt

"""
file_MCN = "./50_75_1000/output_MCN_DAG_50_75_1000.txt"
file = open(file_MCN,"r")

list_prob_time = []
list_time = []

for line_db in file:
    #print(line_db)
    list_line = line_db.split()
    #print(list_line)
    if len(list_line) > 0 and str(list_line[0]) == "Prob":
        proba = str(list_line[3])
        time = str(list_line[7])
        #if proba[0] == "0" or proba[0] == "1": 
        arg = str(list_line[1][4:-1])

        list_prob_time.append([str(arg),float(proba),float(time)])

#print(list_prob_time)
#print(len(list_prob_time))

pickle.dump(list_prob_time,open("MCN_proba_time_50_75_1000", "wb"))
    
file.close()

"""
#print(len(list_prob_time))



#dico_MC = {}

List_diff = []
List_error = []

avg_diff = 0
avg_error = 0

dico_MC = pickle.load(open("MC_MCN_50_75_1000", "rb"))
list_prob_time = pickle.load(open("MCN_proba_time_50_75_1000", "rb"))

cpt = 0

compteur = 0

for i in range(1,1001):
    print("....................")
    print(i)
    print("....................")

    
    #dico_MC[i] = {}

    file_AF = "./50_75_1000/DAG_50_75_"+str(i)+".txt"
    AF = open(file_AF,"r")
    dico_Arg = {}
    dico_Att = {}
    dico_lvl = {}

    for line in AF:
        if line[0:3] == "arg":
            l1 = line.partition("(")
            l2 = l1[2].partition(")")
            dico_Arg[l2[0]] = 1
            dico_lvl[l2[0]] = 0
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
    

    for a in dico_Arg:
        #print(a)

        if cpt > len(list_prob_time) - 1:
            break

        
        arg,p,t = list_prob_time[cpt]
        #print(arg)

        #if a != arg:
        #    print("ERRREUUUUUR !!!!")

        #approximation = MC.approx(str(arg), t, dico_Att, dico_Arg)
        #dico_MC[i][arg] = [approximation,t]

        [approximation,t] = dico_MC[i][arg]

        diff = math.sqrt((p - approximation)**2)

        List_diff.append(diff)
        avg_diff += diff

        #if diff == 0:
        #    List_error.append(float(p)*100)
        #    avg_error += float(p)*100
        #diff/p
        if p == 0:
            List_error.append(100)
            avg_error += 100
        else:
            if float(diff/p) == 1:
                compteur += 1
                #print(diff)
                #print(p)

            List_error.append(float(diff/p)*100)
            avg_error += float(diff/p)*100
        
        cpt += 1

    #pickle.dump(dico_MC, open("MC_MCN_50_75_1000", "wb"))
    #pickle.dump(List_diff, open("diff_MC_MCN_50_75_1000", "wb"))
    #pickle.dump(List_error, open("error_MC_MCN_50_75_1000", "wb"))

#print(cpt)


avg_diff = avg_diff/len(List_diff)    
avg_error = avg_error/len(List_error)       

print(f"avg diff = {avg_diff}")
print(f"avg error = {avg_error}")
print(f"compteur = {compteur}")

#avg diff = 0.05970092489133453
#avg error = 24.268246360975578
#compteur = 3309

#NOOON/ avg diff = 0.05970092489133453
#NOOON/ avg error = 51.14895834810338

#pickle.dump(dico_MC, open("MC_MCN_50_75_1000", "wb"))
#pickle.dump(List_diff, open("diff_MC_MCN_50_75_1000", "wb"))
#pickle.dump(List_error, open("error_MC_MCN_50_75_1000", "wb"))


#avg diff = 0.05970092489133453
#avg error = 51.14895834810338

#dico_MC = pickle.load(open("MC_MCN_50_75_1000", "rb"))
#List_diff = pickle.load(open("diff_MC_MCN_50_75_1000", "rb"))
#List_error = pickle.load(open("error_MC_MCN_50_75_1000", "rb"))

l_error_small = []
for err in List_error:
    if err <= 100 and err>=0:
        l_error_small.append(err)

#pickle.dump(dico_diff, open("Diff_SCN", "wb"))

list_Time = []

for i in range(1,1001):
    for key,val in dico_MC[i].items():
        #print(dico_MC[i])
        #print(key)
        #print(val)
        list_Time.append(val[1])

l_time_small = []
for t in list_Time:
    if t <= 10:
        l_time_small.append(t)



figure = plt.figure(figsize = (10, 10))
plt.gcf().subplots_adjust(left = 0.2, bottom = 0.2,
                       right = 0.7, top = 0.7, 
                       wspace = 0.5, hspace = 0.7)

axes = figure.add_subplot(1, 2, 1)
plt.hist(List_diff,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
#plt.title("Histogram of 1000 MCN (50 arg, 75 att)")
plt.xlabel('Error difference')
plt.yticks([0,1000,3000,5000,10000,15000,20000,25000,30000])
plt.ylabel('Number of occurrence')


axes = figure.add_subplot(1, 2, 2)
#plt.hist(List_error,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
plt.hist(l_error_small,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
#plt.title("Histogram of 1000 MCN (50 arg, 75 att)")
plt.yticks([0,1000,3000,5000,10000,15000,20000])
plt.xlabel('Error percentage')
plt.ylabel('Number of occurrence')


#axes = figure.add_subplot(3, 2, 2)
#plt.boxplot(List_diff)
#plt.ylabel('Error difference')
#plt.title("Boxplot of 1000 graphs - MCN 50,75")
#plt.gca().xaxis.set_ticklabels(['SCN 500'])

#axes = figure.add_subplot(3, 2, 4)
#plt.boxplot(List_error,1)
#plt.ylabel('Error percentage')
#plt.title("Boxplot of 1000 graphs - MCN 50,75")

#axes = figure.add_subplot(3, 2, 5)
#plt.hist(list_Time,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
#plt.hist(l_time_small,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
#plt.title("Histogram of 1000 graphs - MCN 50,75")
#plt.xticks([0,1,2,3,4,5,6,7,8,9,10])
#plt.xlabel('Time (seconds)')
#plt.ylabel('Number of occurrence')

plt.show()

#plt.plot(List_class,List_diff,"-ob") # ob = type de points "o" ronds, "b" bleus
#plt.ylabel('Percentage of error')
#plt.xlabel('Graph class based on number of attacks')
#plt.show()

