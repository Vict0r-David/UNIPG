# Victor DAVID
# Created 03/10/2022

import math 

def grad(dico_Arg,dico_Att):
    new_arg = {}
    for arg,w_ar in dico_Arg.items():
        deg = 1
        for att,w_at in dico_Att.items():
            if arg == att[1]:
               deg = deg * (1-(dico_Arg[att[0]]*(-w_at)))
        new_arg[arg] = deg
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
        print(dico_Arg)
        dico_Arg_pred = dico_Arg
        dico_Arg = grad(dico_Arg, dico_Att)
    return dico_Arg

dico_Arg = {"a":1, "b":1, "c":1, "d":1}
dico_Att = {"ab":-0.6, "bc": -0.3, "db":-0.2, "ad":-0.5, "ca":-0.2}
threshold = 0.01

#dico_Arg = {"a":1, "b":1, "c":1, "d":1, "e":1}
#dico_Att = {"ac":-0.3, "bc": -0.9, "ce":-0.4, "de":-0.3}
#threshold = 0.01

print(semantics(dico_Arg, dico_Att, threshold))