#!/usr/bin/env python3
"""
In this example, we'll start by defining a 3D pyramid as the original model.
The pyramid will have 4 sides and a bottom (or 5 facets/constraints).

We'll then reduce it to a 2D model ignoring the X-direction and see that it
looks like a triangle with three edges and three vertices.

We'll also reduce it to a 2D model ignoring the Y-direction and see that it
looks like a square with four edges and four vertices.

We'll generate plots and graphs for each reduced problem into the
simple_example subfolder
"""
import subprocess
from optlang import *
from fea import flux_envelope_analysis as fea
from fea.plot import plot, generate_graphviz

#### Start by defining our original pyramid ####
model = Model(name='Pyramid')
x,y,z = (Variable('x'),Variable('y'),Variable('z'))
model.add([x,y,z])
model.add(Constraint(y,lb=0,name='base'))
model.add(Constraint(-x+y,ub=1,name='left_wall'))
model.add(Constraint(x+y,ub=1,name='right_wall'))
model.add(Constraint(-z+y,ub=1,name='front_wall'))
model.add(Constraint(z+y,ub=1,name='back_wall'))

#### Let's solve for what it looks like in 2D from the front ####
front_view = fea(model,[x,y])

# Get and save the plot
front_figure, front_plot = plot(front_view, x='x', y='y')
front_figure.savefig("./simple_example/front_view_plot.svg")

# Get graphviz input and generate the image
front_graph = generate_graphviz(front_view)
proc = subprocess.Popen(["dot","-Tpng","-o",'./simple_example/front_view_graph.png'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output = proc.communicate(input=front_graph.encode())[0]


#### Let's solve for what it looks like in 2D from the bottom ####
front_view = fea(model,[x,z])

# Get and save the plot
front_figure, front_plot = plot(front_view, x='x', y='z')
front_figure.savefig("./simple_example/bottom_view_plot.svg")

# Get graphviz input and generate the image
front_graph = generate_graphviz(front_view)
proc = subprocess.Popen(["dot","-Tpng","-o",'./simple_example/bottom_view_plot.png'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output = proc.communicate(input=front_graph.encode())[0]