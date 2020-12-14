#!/usr/bin/env python3
import numpy as np
import matplotlib
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

def _figure_plot(fea,plot=None,figure=None):
    if fea.N>3 or fea.N<2:
        raise ValueError('Cannot plot '+str(fea.N)+' dimensions')
    
    if plot is None:
        if figure is None:
            figure=plt.figure()
        
        labels=[v.name for v in fea._variables]
        if fea.N==3:
            plot=figure.add_subplot(111,projection='3d')
            plot.set_xlabel(labels[0])
            plot.set_ylabel(labels[1])
            plot.set_zlabel(labels[2])
        elif fea.N==2:
            plot=figure.add_subplot(111)
            plot.set_xlabel(labels[0])
            plot.set_ylabel(labels[1])

    return figure, plot

def plot_edges(fea,plot=None,figure=None):
    """
    Plot the edges of an FEA Lattice Graph

    Arguments:
        fea: Lattice Graph/FEA solution to plot
        plot: [Optional] Matplotlib plot
        figure: [Optional] Matplotlib figure
    Returns:
        tuple of Matplotlib figure and plot
    """
    figure, plot=_figure_plot(fea,plot,figure)

    for edge in fea.get_nodes_of_level(1,real=True):
        verts=[tuple(v.point) for v in fea.get_vertices(edge)]
        if len(verts)==2:
            plot.plot(*np.transpose(verts), linestyle='-', linewidth=1, color="k")

    return figure, plot

def plot_vertices(fea,plot=None,figure=None,label=False):
    """
    Plot the vertices of an FEA Lattice Graph

    Arguments:
        fea: Lattice Graph/FEA solution to plot
        plot: [Optional] Matplotlib plot
        figure: [Optional] Matplotlib figure
        label: [Default=False] Whether to label the vertices with their ID
    Returns:
        tuple of Matplotlib figure and plot
    """
    figure, plot=_figure_plot(fea,plot,figure)

    verts=[v for v in fea.get_vertices()]
    
    plot.scatter(*np.transpose([v.point for v in verts]),c="k")

    if label:
        for v in verts:
            if fea.N==3:
                plot.text(x=v.point[0],y=v.point[1],z=v.point[2],s=str(v.id))
            else:
                plot.text(x=v.point[0],y=v.point[1],s=str(v.id))

    return figure, plot

def label_facets(fea,plot=None,figure=None):
    """
    Label the facets of an FEA Lattice Graph.

    Really only useful when combined with plot_edges

    Arguments:
        fea: Lattice Graph/FEA solution to plot
        plot: [Optional] Matplotlib plot
        figure: [Optional] Matplotlib figure
    Returns:
        tuple of Matplotlib figure and plot with the facets labeled
    """
    figure, plot=_figure_plot(fea,plot,figure)

    legends=[]
    points=[]
    for node in fea.get_facets():
        facet=next(iter(node))

        loc=np.mean([v.point for v in fea.get_vertices(node)],axis=0)
        if fea.N==3:
            points.append(plot.scatter(xs=[loc[0]],ys=[loc[1]],zs=[loc[2]]))
        else:
            points.append(plot.scatter(x=[loc[0]],y=[loc[1]]))
        
        legends.append(repr(facet))


    plot.legend(points,legends)

    return figure, plot

__outputs__ = [plot_edges, plot_vertices, label_facets]