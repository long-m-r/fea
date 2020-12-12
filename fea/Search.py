#!/usr/bin/env python3
import logging
log=logging.getLogger('fea.search')

import numpy as np
from .Halfspace import Halfspace
from .util import lstsq
from .VWrapper import VWrapper


class Search:
    """Class for performing a search and calculating the bounding halfspace (if possible)
    """
    _multiplier=10
    _max_iterations=10

    def __init__(self,model,vars,eps=10**-6,clone=True):
        self.interface=model.interface
        if clone:
            self.m=self.interface.Model.clone(model)
        else:
            self.m=model
        self.v=[VWrapper(v,self.m) for v in vars]
        self.eps=eps
        self.n=len(self.v)

        self.H=None
        self.H_cons=None
        self.O=None

    def deactivate(self):
        self.m.remove(self.H_cons)

    def activate(self):
        self.m.add(self.H_cons)
        self.m.objective=self.interface.Objective(np.dot(self.O,self.vexpr))

    def set(self,obj,hps):
        # Deactivate previous constraints
        if self.H_cons is not None:
            self.deactivate()
        
        # Update the search
        self.H=list(hps)
        # self.H_cons=[h._ol_new_constraint(self.vexpr,self.interface,eps=self.eps) for h in self.H]
        self.H_cons=[h._ol_new_constraint(self.m,eps=self.eps) for h in self.H]
        self.O=obj/np.linalg.norm(obj)

        # Go ahead and activate the new problem
        self.activate()

    
    @property
    def vexpr(self):
        return [v.expr for v in self.v]

    @property
    def vp(self):
        return [v.primal for v in self.v]
    @property
    def vd(self):
        return [v.dual for v in self.v]
    
    def Hd(self,i=None):
        if i is not None:
            return self.H[i]._ol_dual(self.m)
        return [h._ol_dual(self.m) for h in self.H]
    
    def Heps(self,i=None):
        if i is not None:
            return self.H[i]._ol_eps(self.m)
        return [h._ol_eps(self.m) for h in self.H]

    def get_solution(self,_i=0):
        self.m.optimize()

        if self.m.status != 'optimal': # Could also check whether =='infeasible'
            log.info('Solver returned status of '+self.m.status)
            if _i<self._max_iterations:
                self.perturb_cons()
                return self.get_solution(_i+1)
            else:
                return False
        log.info('Optimal solution obtained at '+str(self.vp))
        return True

    # Randomly (unless otherwise specified) shift a halfspace constraint
    def perturb_cons(self,index=None):
        if index is None:
            index=np.random.randint(len(self.H))
        
        eps=np.random.uniform(high=self.Heps(index))
        
        if log.getEffectiveLevel()<=10: # Somewhat expensive logging operation
            log.info('Changing '+str(self.H[index])+' to have eps='+str(eps))
        
        rhs=self.H[index].rhs+eps

        self.H_cons[index].ub=None
        self.H_cons[index].lb=rhs
        self.H_cons[index].ub=rhs
        
    # Returns a real Facet or None if it couldn't find anything!
    def bounding_halfspace(self):
        ssemax=(self.eps**2)*self.n

        A=[self.O]
        b=[-1]

        A1_base=np.array([h.norm for h in self.H]+[self.O])
        b1_base=np.array([h._ol_rhs(self.m) for h in self.H]+[self.m.objective.value])

        # If we don't have enough constraints/duals to fully determine the system!
        if len(b1_base)<self.n-1:
            log.info('Insufficient Equations')
            return self.psuedo_halfspace()

        # Loop Through Shadow Prices
        for i,(h,hd) in enumerate(zip(self.H,self.Hd())):
            log.info('Considering '+str(h)+" ~ "+str(hd))
            # TODO: We may not want/need this if statement. A basic halfspace still gives information (i.e. it's perpendicular)

            A1=np.copy(A1_base)
            b1=np.copy(b1_base)

            b1[i]+=self._multiplier
            b1[-1]+=self._multiplier*hd

            try:
                newA=lstsq(A1,b1,self.eps)-self.vp
                A+=[newA]
                b+=[0]
            except ValueError:
                log.info('LstSq Error.')

        # If we can't get enough values to fully determine the system! (This step should only occur when facet duals are 0, which they really shouldn't be)
        if len(A)<self.n:
            log.info('Insufficient Duals')
            return self.psuedo_halfspace()
        # Now, solve for the overall solution
        try:
            log.info('Solving for bounding facet')
            nh=Halfspace(self.vexpr,lstsq(A,b,self.eps),self.vp,eps=self.eps,interface=self.interface)

            return nh

        except ValueError:
            # TODO: Could this be orthoganol?
            log.error('Could not find the bounding facet.')
            return self.psuedo_halfspace()




    # Returns a psuedo facet
    def psuedo_halfspace(self):
        return Halfspace(self.vexpr,-self.O,self.vp,real=False,eps=self.eps,interface=self.interface,required=set(self.H))
