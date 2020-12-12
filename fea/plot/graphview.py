#!/usr/bin/env python3

def _generate_graphview_node(nd):
    # fs="|".join([str(f) for f in nd[0]])
    # if len(fs)>0:
        # fs+=" | "
    # return "\""+repr(nd[0])+"\" ["+"label=\"{"+fs+str(nd[1].get('created','0'))+"}\","+"shape=\"record\","+"color="+"blue"*nd[0].real+"red"*(not nd[0].real)+"];"
    
    if nd[1].get('complete',False):
        return "\""+repr(nd[0])+"\" ["+"label=\"{<in> "+repr(nd[0])+"|"+"|".join([repr(f) for f in nd[0]])+"| <out> "+str(nd[1].get('trace','?'))+"}\","+"shape=\"record\","+"color=blue];"
    elif nd[0].real:
        return "\""+repr(nd[0])+"\" ["+"label=\"{<in> "+repr(nd[0])+"|"+"|".join([repr(f) for f in nd[0]])+"| <out> "+str(nd[1].get('trace','?'))+"}\","+"shape=\"record\","+"color=gray];"
    else:
        return "\""+repr(nd[0])+"\" ["+"label=\"{<in> "+repr(nd[0])+"|"+"|".join([repr(f) for f in nd[0]])+"| <out> "+str(nd[1].get('trace','?'))+"}\","+"shape=\"record\","+"color=red];"

def _generate_graphview_edge(nd):
    searched=nd[2].get('searched',-1)
    if searched>=0:
        return "\""+repr(nd[0])+"\" -> \""+repr(nd[1])+"\" [label=\""+str(searched)+"\", color=\"green\", len=1.00];"
    else:
        return "\""+repr(nd[0])+"\" -> \""+repr(nd[1])+"\";"


def generate_graphview(fea,only_real=False):
    graph=[list() for n in range(fea.N+1)]
    for n in fea.nodes_iter(data=True):
        if n[0].real>=only_real:
            graph[n[0].level].append(_generate_graphview_node(n))

    edges=[]
    for n in fea.edges_iter(data=True):
        if n[0].real>=only_real and n[1].real>=only_real:
            edges.append(_generate_graphview_edge(n))

    retval="digraph G {\n"
    for i,l in enumerate(graph):
        retval+="subgraph cluster"+str(i)+" {\nlabel = \"level "+str(i)+"\";\ncolor = \"gray\";\n"+"\n".join(l)+"\n}\n"
    retval+="\n".join(edges)

    retval+="\n}"

    return retval