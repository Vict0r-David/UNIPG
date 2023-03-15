import networkx as nx
import random

def generateSCN(k,a):
    i = 1
    s = "arg("+ str(a) +"1)."
    liste_arg = [s]
    liste_att = []
    while i < k:
        r = random.randint(1,2)
        for j in range(i+1,i+r+1):
            w = round(random.uniform(0.1,1),2)
            att = "att("+str(a)+str(j)+","+str(a)+str(i)+"):"+str(-w)+"."
            liste_att.append(att)
            liste_arg.append("arg("+str(a)+str(j)+").")
        i += r
    return liste_arg,liste_att

#print(generateSCN(20))

def complex(k,lettres):
    liste_arg = []
    liste_att = []
    for a in lettres:
        l_ar,l_at = generateSCN(k,a)
        liste_arg.append(l_ar)
        liste_att.append(l_at)
    arg = []
    att = []
    arg = arg + liste_arg[0]
    att = att + liste_att[0]
    
    for i in range(1,len(lettres)):
        arg = arg + liste_arg[i]
        att = att + liste_att[i]
        #for j in range(1,k):
            #r = random.randint(2,2)
            #list_r = []
            #for m in range(r):
        r = random.randint(1,k)    
        r2 = random.randint(1,k)
            #while r2 in list_r: 
            #    r2 = random.randint(1,k)
            #list_r.append(r2)
        w = round(random.uniform(0.1,1),2)
        attacks = "att("+str(lettres[i])+str(r)+","+str(lettres[0])+str(r2)+"):"+str(-w)+"."
        att.append(attacks)
    #arg = arg + liste_arg[i+1]
    

    return arg,att

#arg,att = complex(18,["a","b","c","e","f"])
#print(arg)
#print(att)
#print(len(arg))
#print(len(att))


for i in range(1,11):

    liste_arg,liste_att=complex(6,["a","b","c","e","f"])

    file = ".\small_SCN_C\SCN_"+str(i)+".txt"

    fichier = open(file, "w")

    for a in liste_arg:
        fichier.write(str(a)) 
        fichier.write("\n")
    fichier.write("\n")

    for att in liste_att:
        fichier.write(str(att)) 
        fichier.write("\n")
    fichier.write("\n")

    fichier.write("# number of nodes = ")
    fichier.write(str(len(liste_arg)))
    fichier.write("\n")

    fichier.write("# number of edges = ")
    fichier.write(str(len(liste_att)))

    fichier.close()


#for att in liste_e:
 #   print(att)
  #  print(round(DAG[att[0]][att[1]]["weight"],2))
#G.in_edges(node)
#G.out_edges(node) 
