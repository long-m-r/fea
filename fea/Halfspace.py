#!/usr/bin/env python3
import logging
log=logging.getLogger('fea.halfspace')

import numpy as np
import itertools

# TODO: Optimize Static Functions (e.g. key) to not require calculation every time
class Halfspace:
    _id_gen = iter(itertools.count())
    
    # Class for Facet Constraint Storage/Calculation
    def __init__(self,variables,norm,point,real=True,eps=10**-6,interface=None,required=set()):
        self.id=next(Halfspace._id_gen)
        self.eps=eps
        self.real=real
        self.required_halfspaces=required
        
        if len(norm) != len(point):
            raise ValueError('Point and norm do not have the same dimensions')

        self.norm=norm
        self.point=point

        if interface is None:
            interface=variables[0].problem

        self.cons=interface.Constraint(np.dot(self.norm,variables), lb=self.rhs+eps, ub=self.rhs+eps, name=self.name)

    # Calculate point distance and whether it is in the halfspace
    def distance(self,point):
        return np.dot(self.norm,np.array(point)-self.point)

    def contains(self,other,eps=None):
        if eps is None:
            eps=self.eps
        return np.abs(self.distance(other))<=eps
    def __contains__(self,other):
        raise NotImplementedError()
        # return self.contains(other)

    # OptLang Constraint Helper Functions
    def _ol_cons(self,model):
        return model.constraints[self.name]
    def _ol_dual(self,model):
        return self._ol_cons(model).dual
    def _ol_rhs(self,model):
        return self._ol_cons(model).ub
    def _ol_eps(self,model):
        return self._ol_rhs(model)-self.rhs
    def _ol_new_constraint(self,model,eps=None):
        # Cloning is MUCH faster than creating the constraint anew (expression interpretation is the time limiting step!)
        cnew=self.cons.clone(self.cons, model=model)
        if eps is not None:
            cnew.lb=eps+self.rhs
            cnew.ub=eps+self.rhs
        return cnew


# Helpers and Magic Methods
    @property
    def key(self):
        # This is the self contained hashable format for the halfspace. Everything is represented here. Hash and equal will require this
        return tuple([self.real]+list(np.round(np.append(self.norm,self.rhs), self.dec)))

    def __hash__(self):
        return self.key.__hash__()
    
    @property
    def name(self):
        # Originally used hash, had collisions
        return str(self.key).replace(' ','')

    def __eq__(self,other):
        try:
            return self.key == other.key
        except AttributeError:
            return self.__hash__() == other.__hash__()

    def __str__(self):
        return self.__repr__()+': '+'+'.join([str(np.round(self.norm[i],self.dec))+"*(v"+str(i)+"-"+str(np.round(self.point[i],self.dec))+")" for i in range(len(self.norm))])+">"+"PSUEDO"*(not self.real)+"=0+"+str(np.round(self.eps*self.real, self.dec))
    
    def __repr__(self):
        return 'Facet('+str(self.id)+')'

    def __len__(self):
        return len(self.norm)
# Properties for formatting/sanitizing inputs
    @property
    def rhs(self):
        return self._rhs

    @property
    def point(self):
        return self._point
    @point.setter
    def point(self,value):
        self._point=value
        self._rhs=np.dot(self.norm,self._point)

    @property
    def norm(self):
        return self._norm
    @norm.setter
    def norm(self,value):
        self._norm=np.array(value)/np.linalg.norm(value)

    @property
    def eps(self):
        return self._eps
    @eps.setter
    def eps(self,value):
        self._eps=min(value,1)
        self._dec=int(max(0,-np.log10(self._eps)))

    @property
    def dec(self):
        return self._dec
    @dec.setter
    def dec(self,value):
        self._dec=max(0,int(value))
        self._eps=10**(-self._dec)

