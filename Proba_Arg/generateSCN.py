import networkx as nx
import random

def generateSCN(k):
    i = 1
    liste_arg = ["arg(a1)."]
    liste_att = []
    while i < k:
        r = random.randint(1,10)
        for j in range(i+1,i+r+1):
            w = round(random.uniform(0.1,1),2)
            att = "att(a"+str(j)+",a"+str(i)+"):"+str(-w)+"."
            liste_att.append(att)
            liste_arg.append("arg(a"+str(j)+").")
        i += r
    return liste_arg,liste_att

#print(generateSCN(20))


for i in range(1,21):

    liste_arg,liste_att=generateSCN(800)

    file = ".\SCN_800\SCN_800_"+str(i)+".txt"

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
