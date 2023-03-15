import math 

file_Approx = ".\DAG_20+\output_Approx_DAG+_20.txt"
file_A = open(file_Approx,"r")

list_approx = []

for line_db in file_A:
    #print(line_db)
    list_line = line_db.split()
    #print(list_line)
    if len(list_line) > 0 and str(list_line[0]) == "Prob":
        proba = str(list_line[3])
        #if proba[0] == "0" or proba[0] == "1": 
        list_approx.append(proba)
    
file_A.close()


file_MCN = ".\DAG_20+\output_MCN_DAG+_20.txt"
file_M = open(file_MCN,"r")

list_mcn = []

for line_db in file_M:
    list_line = line_db.split()
    if len(list_line) > 0 and str(list_line[0]) == "Prob":
        proba = str(list_line[3])
        #if proba[0] == "0" or proba[0] == "1": 
        list_mcn.append(proba)

file_M.close()

avg_dist = 0
print(len(list_approx))
print(len(list_mcn))
for i in range(0,len(list_approx)):
    avg_dist += abs(float(list_approx[i]) - float(list_mcn[i]))
print(avg_dist/len(list_approx))