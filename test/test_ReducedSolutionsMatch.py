#!/usr/bin/env python3
import unittest
import numpy as np
from optlang import *
from fea import flux_envelope_analysis as fea
from itertools import permutations


class ReducedSolutionsMatch(unittest.TestCase):
    def setUp(self):
        # Parameters for our full/original model
        self.original_dims=6
        self.original_cons=2*self.original_dims
        self.var_limit=10
        self.nsteps=10

        # Construct the original model
        print('Constructing random '+str(self.original_dims)+'D model')
        self.model=Model(name='Random')

        ## Add variables
        self.original_vars=[]
        for i in range(self.original_dims):
            self.original_vars.append(Variable(chr(65+i),ub=self.var_limit,lb=-self.var_limit))
        self.model.add(self.original_vars)

        ## Add random constraints
        for i in range(self.original_cons):
            val=np.random.rand(self.original_dims)
            val=val/np.linalg.norm(val)
            cns=(np.random.random()-.5)*2*self.var_limit
            if cns>0:
              self.model.add(Constraint(np.dot(self.original_vars,val),ub=cns,name='C'+str(i)))
            else:
              self.model.add(Constraint(np.dot(self.original_vars,val),lb=cns,name='C'+str(i)))

        # Create the reduced solution via FEA
        print('Solving for reduced '+str(self.original_dims//2)+'D model')
        self.subset=self.original_vars[0:self.original_dims//2]
        self.reduced=fea(self.model,self.subset)
        self.reduced_optlang=self.reduced.to_optlang_model(False)
        self.facets=[next(iter(a)) for a in self.reduced.get_facets()]

    def stepwise_facet_check(self, model,variables,vals=[]):
        """recursive function for checking facets match the model.

        Will repeat for nsteps increments per dimension
        """
        cv=variables[0]
        self.model.objective=Objective(cv, direction='max')
        self.model.optimize()
        vmax=cv.primal
        self.model.objective.direction='min'
        self.model.optimize()
        vmin=cv.primal

        if len(variables)==1:
            # We're at the last point, so check the distance
            md=min([f.distance(vals+[vmin]) for f in self.facets]+[f.distance(vals+[vmax]) for f in self.facets])

            # Ensure the distance is within the margin of error
            self.assertTrue(md > -self.reduced.EPS)
            self.assertTrue(md < self.reduced.EPS)
        else:
            # We have more dimensions to reduce, recurse
            vlb=cv.lb
            vub=cv.ub
            for v in np.linspace(vmin,vmax,self.nsteps):
                cv.set_bounds(v,v)
                self.stepwise_facet_check(self.model,variables[1:],vals=vals+[v])
            cv.set_bounds(vlb,vub)

    def test_randomRaysHaveSameSolution(self):
        """Ensure FEA solution matches random original problem by finding the solution to random rays"""
        random_trys=100
        print('Checking '+str(random_trys)+' random rays inside the original and reduced model to ensure they match')

        for r in range(random_trys):
            # Get a random linear objective function for this using only our subset variables
            lc={k:v for k,v in zip(self.subset,np.random.rand(self.original_dims//2)-.5)}

            # Solve the original model for that objective and get values
            self.model.objective.set_linear_coefficients(lc)
            self.reduced_optlang.objective.set_linear_coefficients(lc)
            self.model.optimize()
            mv=self.model.objective.value
            ptm=np.array([self.model.primal_values[v.name] for v in self.subset])

            # Solve the reduced model for that objective and get values
            self.reduced_optlang.optimize()
            rv=self.reduced_optlang.objective.value
            ptr=np.array([self.reduced_optlang.primal_values[v.name] for v in self.subset])

            # We should get the same point
            self.assertTrue(np.linalg.norm(ptm-ptr) <= self.reduced.EPS**2)
            # We should get the same value
            self.assertTrue(np.abs(mv-rv) <= self.reduced.EPS)

    def test_facetsMatchOriginalModel(self):
        """Ensure that facets from reduced model are contained within the original model.

        This test is recursive so use a helper method
        """
        print('Checking that reduced facets match the original model via stepwise approximation')
        self.stepwise_facet_check(self.model, self.subset)

    def test_furtherReductionsMatchFromBothModels(self):
        """Reduce the original model and our reduced model to a 2-dimensional model. These should match"""
        print('Checking that all possible 2D models created from the original and reduced models match')
        
        for combo in permutations(self.subset,2):
          # Solve both for the current 2D problem
          a=fea(self.model,combo)
          b=fea(self.reduced_optlang,combo)
          
          # They should have the same facets
          disjoint_facets=a.get_facets()^b.get_facets()
          self.assertEqual(len(disjoint_facets),0)
          
          # They should have the same vertices
          disjoint_vertices=a.get_vertices()^b.get_vertices()
          self.assertEqual(len(disjoint_vertices),0)
