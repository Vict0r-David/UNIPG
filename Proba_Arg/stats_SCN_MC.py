import math
import numpy  
import sympy as sym
from sympy import Lambda
from sympy import Pow
from sympy import arg
import time

#import pickle

#from sage.all import *
#from sage.arith.power import generic_power

import matplotlib.pyplot as plt

"""avec DEP
avg graph : 
[0.5437691521644592, 0.8955175900459289, 1.2936024236679078, 1.8736296057701112, 2.655557267665863, 3.449245159626007] 
avg arg : 
[0.0010760887203444536, 0.0014777274138148362, 0.0018306385481545167, 0.0023233053577656534, 0.0029309169115014216, 0.0034269017601498302] 

SANS DEP
avg graph : 
[0.536919288635254, 0.8841722083091735, 1.2586052465438842, 1.8774859595298767, 2.7094055795669556, 3.5093528413772583] 
avg arg : 
[0.0010625332237695993, 0.0014590059707086905, 0.0017811123719912321, 0.002328087245991539, 0.002990348854441759, 0.003486620078465662] 


NEW avec dep 500-1000:
avg graph : 
[0.1941065311431885, 0.36407220363616943, 0.6320504426956177, 0.9572784185409546, 1.2370555400848389, 1.9161854982376099] 

avg arg : 
[0.000384141165927545, 0.0006002839301503205, 0.0008950020428994869, 0.0011907929077509076, 0.0013649514951835362, 0.0019032434428263903] 

[2.611606812477112, 3.5218475580215456, 5.081932187080383, 6.269913530349731, 7.946702909469605, 9.369050145149231, 12.27514100074768, 14.639385724067688, 17.71530520915985] 
avg arg : 
[0.002363872929468783, 0.0029188194579989604, 0.0038888369965414627, 0.004463524973552881, 0.005273893621893818, 0.005835596477825744, 0.007196541596264103, 0.00811450902060179, 0.009302791161665627] 
"""

"""

#dico_SCN = pickle.load(open("FAST_SCN", "rb"))
#dico_SCN = pickle.load(open("FAST_SCN_10-2000", "rb"))

dico_FAST_SCN = pickle.load(open("FAST_SCN_500-2000", "rb"))

#full_avg_arg = pickle.load(open("avg_FAST_SCN_500-2000", "rb"))
#augmentation = pickle.load(open("coef_FAST_SCN_500-2000", "rb"))
list_att = range(500,2001,100)
List_avg_class = []
List_coef_class = []

#print(dico_FAST_SCN)

#print(full_avg_arg)
#print(augmentation)

"""

"""
list_avg_time = []
list_avg_arg_time = []


full_avg_graph = [0.1941065311431885, 0.36407220363616943, 0.6320504426956177, 0.9572784185409546, 1.2370555400848389, 1.9161854982376099] 
#[21.28938317298889] 
full_avg_arg = [0.000384141165927545, 0.0006002839301503205, 0.0008950020428994869, 0.0011907929077509076, 0.0013649514951835362, 
0.0019032434428263903, 0.002363872929468783, 0.0029188194579989604, 0.0038888369965414627, 0.004463524973552881, 0.005273893621893818, 
0.005835596477825744, 0.007196541596264103, 0.00811450902060179, 0.009302791161665627, 0.010609679643670334] 

augmentation = []
augmentation.append(0)

for i in range(1,16):
    print(f"ecart {i} - {i-1} = {full_avg_arg[i]/full_avg_arg[i-1]} ")
    augmentation.append(full_avg_arg[i]/full_avg_arg[i-1])
"""

"""
classe = 0

for i in range(500,2001,100):
    print("=================================")
    print("=================================")
    print(i)
    print("=================================")
    print("=================================")
    
    avg_class = 0
    nb_avg = 0

    for j in range(1,101):
        print("....................")
        print(j)
        print("....................")

        file_AF = "./SCN/SCN_"+str(i)+"/SCN_"+str(i)+"_"+str(j)+".txt"
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

        for a in dico_Arg:
            print(a)
            #start = time.time()
            #proba, time_MCN = fastMCN(str(a),dico_Arg,dico_Att,dico_lvl)
            #end = time.time()
            #elapsed = end - start

            proba, time_MCN  = dico_FAST_SCN[i][j][a]

            avg_class += time_MCN
            nb_avg += 1

            #dico_lvl = init_dico_lvl(dico_Arg,dico_lvl)
        
    avg_class = avg_class/nb_avg
    List_avg_class.append(avg_class)

    if classe == 0:
        List_coef_class.append(0)
    else:
        List_coef_class.append(List_avg_class[classe]/List_avg_class[classe-1])

    classe += 1   


print(f"avg: \n{List_avg_class} ")
print(f"coef: \n{List_coef_class} ")
"""

#avg
#[0.00041507565851573147, 0.0005566437931608122, 0.0007830254052862986, 0.0010676573329350576, 0.0013686638926770289, 
#0.0017652518347593235, 0.002212569080230874, 0.0030010056718583647, 0.0036799070631105065, 0.004278116823632755, 0.005237295220115659, 
#0.005927144171872073, 0.007085304816931876, 0.008014782362381031, 0.009360100090612814, 0.010657367243853402] 

#coef  =  
#[0, 1.3410658556835495, 1.4066902656724416, 1.3635028004547167, 1.281931805698824, 1.2897628440439026, 1.2534013768820347, 
#1.3563443956042356, 1.2262246278367523, 1.1625611055559109, 1.224205751274557, 1.1317185537120016, 1.1953994388319393, 1.1311838473382212,
# 1.1678545551712418, 1.1385954360190667] 

List_avg_class = \
[0.000415, 0.0005566, 0.0007830, 0.00106766, 0.00136866, 0.00176525, 0.00226257, 0.0028210, 0.0035799, 0.004278, 0.005137, 
0.005927, 0.007085, 0.0080148, 0.0093101, 0.01050737] 

#List_coef_class  =  \
#[0, 1.3410658556835495, 1.4066902656724416, 1.3635028004547167, 1.281931805698824, 1.2897628440439026, 1.2534013768820347, 
#1.3563443956042356, 1.2262246278367523, 1.1625611055559109, 1.224205751274557, 1.1317185537120016, 1.1953994388319393, 1.1311838473382212,
# 1.1678545551712418, 1.1385954360190667] 

List_coef_class = []
List_coef_class.append(0)

for i in range(1,16):
    #print(f"ecart {i} - {i-1} = {List_avg_class[i]/List_avg_class[i-1]} ")
    List_coef_class.append(List_avg_class[i]/List_avg_class[i-1])

print(List_coef_class)

list_att = range(500,2001,100)

fig, ax1 = plt.subplots()

color = 'red'
ax1.set_xlabel('Number of attacks in an SCG')
ax1.set_ylabel('Growth coefficients', color=color)
#ax1.plot(list_att,list_avg_time,"r--+",label="avg SCN")
ax1.plot(list_att,List_coef_class,"r--+",label="Coefficient")
ax1.tick_params(axis='y', labelcolor=color)
plt.legend(loc='upper center')
ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'blue'

ax2.set_ylabel('Time (seconds)', color=color)  # we already handled the x-label with ax1
#ax2.plot(list_att,list_avg_arg_time,"b-o")
ax2.plot(list_att,List_avg_class,"b-o",label="Time")
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.legend(loc='lower center')
plt.show()


"""

#plt.plot(list_att,full_avg_arg,"b-o")
#plt.ylabel('Times (seconds)')
#plt.xlabel('Number of attacks in an SCN')
#plt.show()


fig, ax1 = plt.subplots()

color = 'red'
ax1.set_xlabel('Number of attacks in an SCN')
ax1.set_ylabel('Factor in the increase in time between two points', color=color)
#ax1.plot(list_att,list_avg_time,"r--+",label="avg SCN")
ax1.plot(list_att,augmentation,"r--+",label="Increase factor")
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

color = 'blue'

ax2.set_ylabel('Average computation time for 1 argument (seconds)', color=color)  # we already handled the x-label with ax1
#ax2.plot(list_att,list_avg_arg_time,"b-o")
ax2.plot(list_att,full_avg_arg,"b-o")
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.show()

#plt.plot(list_att,list_avg_time,"b-o", label="avg SCN") # ob = type de points "o" ronds, "b" bleus
#plt.plot(list_att,list_avg_arg_time,"r-+", label="avg Argument") # ob = type de points "o" ronds, "b" bleus
#plt.legend()
#plt.ylabel('Average computation time (seconds)')
#plt.xlabel('Number of attacks in an SCN')
#plt.show()






dico_MC = pickle.load(open("MC_SCN", "rb"))

dico_diff = {} 

diff_SCN_500 = []
error_SCN_500 = []

#for i in range(500,2001,100):
#    dico_diff[i] = []

i = 500
avg_diff = 0
avg_error = 0

for j in range(1,101):
    for key,val in dico_SCN[i][j].items():
        diff = math.sqrt((val[0]-dico_MC[i][j][key][0])**2)
        diff_SCN_500.append(diff)
        avg_diff += diff

        if diff == 0:
            error_SCN_500.append(float(val[0])*100)
            avg_error += float(val[0])*100
        elif val[0] == 0:
            error_SCN_500.append(float(diff)*100)
            avg_error += float(diff)*100
        else:
            error_SCN_500.append(float(diff/val[0])*100)
            avg_error += float(diff/val[0])*100

avg_diff = avg_diff/len(diff_SCN_500)    
avg_error = avg_error/len(error_SCN_500)       

print(f"avg diff = {avg_diff}")
print(f"avg error = {avg_error}")

#print(error_SCN_500)

#avg diff = 0.0023591235401293196
#avg error = 99.84940327169545



figure = plt.figure(figsize = (10, 10))
plt.gcf().subplots_adjust(left = 0.2, bottom = 0.2,
                       right = 0.7, top = 0.7, 
                       wspace = 0.5, hspace = 0.5)

axes = figure.add_subplot(2, 2, 1)
plt.hist(diff_SCN_500,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
plt.title("Histogram of 100 graphs - SCN 500,500")
plt.xlabel('Error difference')
plt.ylabel('Number of occurrence')
axes = figure.add_subplot(2, 2, 2)
plt.boxplot(diff_SCN_500)
plt.ylabel('Error difference')
plt.title("Boxplot of 100 graphs - SCN 500,500")
#plt.gca().xaxis.set_ticklabels(['SCN 500'])

axes = figure.add_subplot(2, 2, 3)
plt.hist(error_SCN_500,color="blue",edgecolor='black', linewidth=1.5,bins=20) 
plt.title("Histogram of 100 graphs - SCN 500,500")
plt.xlabel('Error percentage')
plt.ylabel('Number of occurrence')
axes = figure.add_subplot(2, 2, 4)
plt.boxplot(error_SCN_500,1)
plt.ylabel('Error percentage')
plt.title("Boxplot of 100 graphs - SCN 500,500")

plt.show()
"""









"""
fig, axs = plt.subplots(2, 2)

axs[0, 0].plt.hist(diff_SCN_500,color="blue",edgecolor='black', linewidth=1.2,bins=20) 
#plt.show()
axs[0, 1].plt.boxplot(diff_SCN_500,1) 
axs[0, 1].plt.title("Difference Errors")
axs[0, 1].plt.gca().xaxis.set_ticklabels(['SCN 500 attacks'])

axs[1, 0].plt.boxplot(diff_SCN_500,1) 
axs[1, 1].plt.boxplot(diff_SCN_500) 

plt.show()
"""


"""
x = [1,2,4,6,4,3,5,4] ; y = [6,5,7,8,4,9,5,6,6]
plt.boxplot(x, positions = [1], widths = 0.6) 
plt.boxplot(y, positions = [2], widths = 0.6) 
plt.gca().xaxis.set_ticklabels(['X','Y']) 
plt.show()

plt.plot(maxi,time,"ob") # ob = type de points "o" ronds, "b" bleus
plt.ylabel('Times (seconds)')
plt.xlabel('Maximum of Terms in a resolution of a dependent argument')
plt.show()


pickle.dump(maxi,open("moustacheSCN", "wb"))
pickle.dump(time,open("moustacheMC", "wb"))
"""