#!/usr/bin/env python3

import numpy as np
import matplotlib
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Line3DCollection

def plot(graph, **kwargs):
    """
    Plots a 2d or 3d graph given an FEA LatticeGraph object (the output of the fea function)

    Arguments:
        graph: FEA Lattice Graph to Plot

    Named Arguments:
        figure: The figure to plot on
        cmap: color map (defaults to inferno)
        x: specify variable name for X-axis
        y: specify variable name for Y-axis
        z: specify variable name for Z-axis
        c: specify variable name for color-axis
        vertices: Default True. Whether or not to plot vertices
        edges: Default True. Whether or not to plot edges

    Returns:
        Tuple of matplotlib figure and plot

    Example:
        solution = fea(model, [x,y,z])
        solfigure, solplot = plot(solution, x='x', y='y', z='z')
        solfigure.savefix('solution_figure.svg')
    """
    figure=kwargs.get('figure',plt.figure())
    cm=kwargs.get('cmap',matplotlib.cm.get_cmap('inferno'))
    
    if kwargs.get('z',None) is not None:
        plot=kwargs.get('plot',figure.add_subplot(111,projection='3d'))
        plot.set_xlabel(str(kwargs.get('x','')))
        plot.set_ylabel(str(kwargs.get('y','')))
        plot.set_zlabel(str(kwargs.get('z','')))
        
        dname=[str(kwargs['x']),str(kwargs['y']),str(kwargs['z'])]
    else:
        plot=kwargs.get('plot',figure.add_subplot(111))
        plot.set_xlabel(str(kwargs.get('x','')))
        plot.set_ylabel(str(kwargs.get('y','')))
        
        dname=[str(kwargs.get('x','')),str(kwargs.get('y',''))]
        
    dmat=np.zeros([len(dname),graph.N])
    vnames=[str(v) for v in graph._variables]
    
    if kwargs.get('c',None) in vnames:
        c=None
        cmat=np.zeros(graph.N)
        cmat[vnames.index(kwargs['c'])]=1
    else:
        c=kwargs.get('c','k')
        cmat=None
        
    for i,n in enumerate(dname):
        dmat[i,vnames.index(str(n))]=1
    
        
    # Plot vertices (if desired)
    if kwargs.get('vertices',True):
        verts=list(graph.get_vertices(real=kwargs.get('real',True),complete=kwargs.get('complete',True)))
        
        if c is None:
            ccurr=[np.dot(cmat,v.point) for v in verts]
        else:
            ccurr=c
        
        sc=plot.scatter(*np.transpose([np.dot(dmat,v.point) for v in verts]),c=ccurr,cmap=cm)
    
       
    # Plot Edges (if desired)
    if kwargs.get('edges',True):
        edges=list(graph.get_nodes_of_level(1,real=kwargs.get('real',True),complete=kwargs.get('complete',True)))
        
        edge_collection=[]
        if c is None:
            ccurr=[]
        else:
            ccurr=c
        
        for e in edges:
            verts=list(graph.get_vertices(e,real=kwargs.get('real',True),complete=kwargs.get('complete',True)))
            if len(verts)==2:
                if c is None:
                    dat_range=np.array([np.dot(dmat,v.point) for v in verts])
                    c_range=np.array([np.dot(cmat,v.point) for v in verts])
                    
                    steps=50
                    
                    
                    points=np.array([np.linspace(r[0],r[1],steps) for r in dat_range.T]).T.reshape(-1,1,len(dname))
                    segs=np.concatenate([points[:-1],points[1:]], axis=1)
                    edge_collection.extend(list(segs))
                    
                    ccurr.extend(np.linspace(c_range[0],c_range[1],steps-1))

                else:
                    edge_collection.append([np.dot(dmat,v.point) for v in verts])


        ax = plot.axes    
        if c is None:
            if len(dname)==3:
                lc=Line3DCollection(edge_collection, array=np.array(ccurr),cmap=cm)
                ax.add_collection3d(lc)
            else:
                lc=LineCollection(np.array(edge_collection), array=np.array(ccurr),cmap=cm)
                ax.add_collection(lc)    
        else:
            if len(dname)==3:
                lc=Line3DCollection(np.array(edge_collection))
                ax.add_collection3d(lc)
            else:
                lc=LineCollection(np.array(edge_collection))
                ax.add_collection(lc)    
            
    return figure,plot

# def fix_clims(f,p,bounds=None,mul=[1,1]):
#     f.show()
#     if bounds is None:
#         print([c.get_clim() for c in p.collections])
#         lims=np.array([c.get_clim() for c in p.collections])
#         bounds=[min((i for i in lims.T[0] if i is not None)),max((i for i in lims.T[1] if i is not None))]
#     bounds=np.multiply(bounds,mul)

#     for c in p.collections:
#         c.set_clim(*bounds)
    
#     f.colorbar(p.collections[0])
#     return bounds

__exports__ = [plot]