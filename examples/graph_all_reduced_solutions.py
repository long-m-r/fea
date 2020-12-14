#!/usr/bin/env python3
"""
This example shows how to generate basic plots and manipulate them with FEA

In this case, we'll generate a 4D model and ruse FEA to find all possible
2D and 3D reduced models.

We'll then solve for the solution spaces by breaking up each solution and
finding the maximum and minimum values in a 10 step grid in each dimension.

Finally we'll generate a graph of the reduced solution showing the nodes
we discovered

Note: If you are using GLPK as your solver your results may vary. I highly
recommend using CPLEX or Gurobi if you are licensed to do so.
"""
import subprocess
import numpy as np
from optlang import *
from itertools import combinations, chain
from fea import flux_envelope_analysis as fea
from fea.plot import generate_graphviz

#### EXAMPLE PARAMETERS ####
# how many dimensions and constraints we want in our original random model
original_dims=4
original_cons=original_dims*2

# upper and lower bounds for our variable
var_limit=10

#### MODEL FORMULATION AND HELPER FUNCTIONS ####

# Generate a Random Model
np.random.seed(1)
model=Model(name='Random')
original_vars=[]
for i in range(original_dims):
    original_vars.append(Variable(chr(65+i),ub=var_limit,lb=-var_limit))
model.add(original_vars)
for i in range(original_cons):
    val=np.random.rand(original_dims)
    val=val/np.linalg.norm(val)
    cns=(np.random.random()-.5)*var_limit
    if cns>0:
      model.add(Constraint(np.dot(original_vars,val),ub=cns,name='C'+str(i)))
    else:
      model.add(Constraint(np.dot(original_vars,val),lb=cns,name='C'+str(i)))

#### SOLVE AND GRAPH ####
# Use all possible 2d and 3d combos
for combo in chain(combinations(original_vars, 2),combinations(original_vars, 3)):
  # Solve FEA model (if using GLPK, expect the 3d ones to be very noisy)
  reduced = fea(model, combo)

  # Get graphviz format
  reduced_graph = generate_graphviz(reduced)

  # Generate the image from the input
  proc = subprocess.Popen(["dot","-Tpng","-o",'./graph_all_reduced_solutions/'+''.join([v.name for v in combo])+'.png'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output = proc.communicate(input=reduced_graph.encode())[0]

