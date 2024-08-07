import networkx as nx
import random


for i in range(1,1501):
        
    G=nx.gnp_random_graph(50,0.1,directed=True)

    DAG = nx.DiGraph([(u,v,{'weight':random.uniform(0.1,1)}) for (u,v) in G.edges() if u<v])

    #DAG.edges =
    #print(DAG.edges)
    l0 = list(DAG.edges)
    #print(type(l0))
    j = 0
    l = []
    while j < 70:
        l.append(l0[j])
        j +=1
    li = list(set(l0) - set(l))

    DAG.remove_edges_from(li)
    #print(f"edge = {DAG.edges}")
    #print(f"l = {l}")
    #DAG.add_edges_from(l)
    #print(f"l = {l}")
    #print(f"edge = {DAG.edges}")

    #print(DAG)
    #print(nx.is_directed_acyclic_graph(DAG))

    liste_n = DAG.nodes
    liste_arg = []
    for n in liste_n:
        liste_arg.append("arg(a"+str(n)+").")
    #print(liste_arg)

    liste_e = DAG.edges
    #dico_att = {}
    liste_att = []
    for att in liste_e:
        liste_att.append("att(a"+str(att[0])+",a"+str(att[1])+"):"+str(-round(DAG[att[0]][att[1]]["weight"],2))+".")
        #dico_att["a"+str(att[0])+"->a"+str(att[1])] = ["a"+str(att[0]),"a"+str(att[1]),-round(DAG[att[0]][att[1]]["weight"],2)]
    #print(liste_att)


    file = "./../../50_70_1500/DAG_50_70_"+str(i)+".txt"

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