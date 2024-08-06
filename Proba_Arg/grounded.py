def list_dico_Att_in(dico_Arg, dico_Att):
    ldico_att_in = {}
    for id_arg in dico_Arg:
        ldico_att_in[id_arg] = []
    for id,att in dico_Att.items():
        ldico_att_in[att[1]] += [id]
    return ldico_att_in


def dlist_Att_in_to_Arg(dico_list_att):
    dico = {}
    for arg,l in dico_list_att.items():
        dico[arg] = []
        for att in l:
            l2 = att.split("-")
            dico[arg].append(l2[0])
    return dico


def list_dico_Att_out(dico_Arg, dico_Att):
    ldico_att_out = {}
    for id_arg in dico_Arg:
        ldico_att_out[id_arg] = []
    for id,att in dico_Att.items():
        ldico_att_out[att[0]] += [id]
    return ldico_att_out

def dlist_Att_out_to_Arg(dico_list_att):
    dico = {}
    for arg,l in dico_list_att.items():
        dico[arg] = []
        for att in l:
            l2 = att.split("->")
            dico[arg].append(l2[1])
    return dico

def grounded(dico_Arg,dico_Att):  
    d_in = list_dico_Att_in(dico_Arg, dico_Att)
    att_in = dlist_Att_in_to_Arg(d_in)
    d_out = list_dico_Att_out(dico_Arg, dico_Att)
    att_out = dlist_Att_out_to_Arg(d_out)
    to_be_in = []
    label = {}
    und_pre = {}
    for x in dico_Arg:
        label[x] = "und"
        und_pre[x] = len(att_in[x]) 
        if und_pre[x] == 0:
            to_be_in.append(x)
 
    while len(to_be_in) > 0:
        x = to_be_in[0]
        to_be_in.remove(x)
        label[x] = "in"
        for y in att_out[x]:
            if label[y] != "out":
                label[y] = "out"
                for z in att_out[y]:
                    if label[z] == "und":
                        und_pre[z] = und_pre[z] - 1
                        if und_pre[z] == 0:
                            to_be_in.append(z)
    output = [] 
    for arg in label:
        if label[arg] == "in":
            output.append(arg)
    return output