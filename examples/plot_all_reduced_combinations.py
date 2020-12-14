#!/usr/bin/env python3
"""
This example shows how to generate basic plots and manipulate them with FEA

In this case, we'll generate a 4D model and ruse FEA to find all possible
2D and 3D reduced models.

We'll then solve for the solution spaces by breaking up each solution and
finding the maximum and minimum values in a 10 step grid in each dimension.

Finally we'll generate a plot of the reduced solution and add lines showing
the results from our grid search (which should line up perfectly).

Note: If you are using GLPK as your solver your results may vary. I highly
recommend using CPLEX or Gurobi if you are licensed to do so.
"""
import numpy as np
import matplotlib
from optlang import *
from itertools import combinations, chain

from fea import flux_envelope_analysis as fea
from fea.plot import plot

#### EXAMPLE PARAMETERS ####
# how many dimensions and constraints we want in our original random model
original_dims=4
original_cons=original_dims*2

# upper and lower bounds for our variable
var_limit=10

# number of steps to use for our stepwise comparison
nsteps=10

# Figure Defaults
matplotlib.rc('font', size=22)

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


def stepwise_check(model,variables):
  """Function for doing a stepwise max/min to find the general shape of the solution space"""
  cv=variables[0]
  model.objective=Objective(cv, direction='max')
  model.optimize()
  vmax=cv.primal
  model.objective.direction='min'
  model.optimize()
  vmin=cv.primal

  if len(variables)==1:
    return [[[vmin],[vmax]]]
  else:
    retval=[]
    vlb=cv.lb
    vub=cv.ub
    for v in np.linspace(vmin,vmax,nsteps):
      cv.set_bounds(v,v)
      retval.extend([[[v]+p for p in s] for s in stepwise_check(model,variables[1:])])
    cv.set_bounds(vlb,vub)
    return retval

#### SOLVE AND PLOT ####
# Use all possible 2d and 3d combos
for combo in chain(combinations(original_vars, 2),combinations(original_vars, 3)):
  # Solve FEA model (if using GLPK, expect the 3d ones to be very noisy)
  reduced = fea(model, combo)

  # Solve stepwise (for comparison)
  stepwise = stepwise_check(model, combo)

  # Plot the reduced solution
  currfigure, currplot = plot(reduced, x=combo[0].name, y=combo[1].name)

  # Add the stepwise stuff
  for step in stepwise:
    currplot.plot(*np.transpose(step),"X:",markersize=5, linewidth=1, color="gray")

  # Save the figure
  currfigure.savefig("./plot_all_reduced_combinations/"+''.join([v.name for v in combo])+".svg")

