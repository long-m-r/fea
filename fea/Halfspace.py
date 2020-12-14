#!/usr/bin/env python3
import logging
log=logging.getLogger('fea.halfspace')

import numpy as np
import itertools

class Halfspace:
    _id_gen = iter(itertools.count())
    
    # Class for Facet Constraint Storage/Calculation
    def __init__(self,variables,norm,point,real=True,eps=10**-6,interface=None,required=set()):
        """
        Class for a facet/halfspace constraint in the LatticeGraph as defined in :class:`fea.LatticeGraph`
        Combinations of Halfspace objects define a :class:`fea.Node`

        They are constructed of a series of variables, a normal unit vector, and a point

        Parameters
        ----------
        variables: array of :class:`optlang.Variable`
            The variables this halfspace is defined over
        norm:   :class:`numpy.ndarray`
            Array indicating the normal vector for this halfspace
        point:  :class:`numpy.ndarray`
            Array indicating a point in space that lies on this halfspace
        real:   bool
            Whether or not this is a real or psuedo-halfspace (one created by the solver for processing)
        eps:    float
            The detection limit/margin of error for the solver
        interface:
            Optlang interface to use. It'll steal it from the variables if not defined
        required: array of :class:`fea.Halfspace`
            Halfspaces which are required for this one to be utilized. Only necessary when real=False
        """
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

    def distance(self,point):
        """
        Calculate the distance between this halfspace and a point

        Parameters
        -----------
        point: :class:`numpy.ndarray`
            A point to check the distance against
        """
        return np.dot(self.norm,np.array(point)-self.point)

    def contains(self,other,eps=None):
        """
        Check whether this halfspace contains a point

        Respects the detection margin, eps

        Parameters
        -----------
        other: :class:`numpy.ndarray`
            A point to check whether it is contained
        eps: float
            The detection limit. Defaults to class value, but can be overriden
        """
        if eps is None:
            eps=self.eps
        return np.abs(self.distance(other))<=eps

    def __contains__(self,other):
        """Should indicate whether this halfspace contains a Node or HalfSpace.

        Not implemented yet"""
        raise NotImplementedError()
        # return self.contains(other)

    # OptLang Constraint Helper Functions
    def _ol_cons(self,model):
        """Private. Get OptLang Constraint for this halfspace"""
        return model.constraints[self.name]
    def _ol_dual(self,model):
        """Private. Get OptLang dual for this halfspace"""
        return self._ol_cons(model).dual
    def _ol_rhs(self,model):
        """Private. Get OptLang upper bound for this halfspace"""
        return self._ol_cons(model).ub
    def _ol_eps(self,model):
        """Private. Get OptLang eps value for this halfspace"""
        return self._ol_rhs(model)-self.rhs
    def _ol_new_constraint(self,model,eps=None):
        """Private. Generate OptLang Constraint for this halfspace"""
        # Cloning is MUCH faster than creating the constraint anew (expression interpretation is the time limiting step!)
        cnew=self.cons.clone(self.cons, model=model)
        if eps is not None:
            cnew.lb=eps+self.rhs
            cnew.ub=eps+self.rhs
        return cnew

# Helpers and Magic Methods
    @property
    def key(self):
        """A string representation of this hashtable. Guaranteed to match for identical halfspaces (within our detection limit, eps)"""
        # This is the self contained hashable format for the halfspace. Everything is represented here. Hash and equal will require this
        return tuple([self.real]+list(np.round(np.append(self.norm,self.rhs), self.dec)))

    def __hash__(self):
        """Use the key to generate a hash"""
        return self.key.__hash__()

    @property
    def name(self):
        """Use this as a distinct name"""
        # Originally used hash, had collisions
        return str(self.key).replace(' ','')

    def __eq__(self,other):
        """Check whether this halfspace and another are identical"""
        try:
            return self.key == other.key
        except AttributeError:
            return self.__hash__() == other.__hash__()

    def __str__(self):
        return self.__repr__()+': '+'+'.join([str(np.round(self.norm[i],self.dec))+"*(v"+str(i)+"-"+str(np.round(self.point[i],self.dec))+")" for i in range(len(self.norm))])+">"+"PSUEDO"*(not self.real)+"=0+"+str(np.round(self.eps*self.real, self.dec))

    def __repr__(self):
        return 'Facet('+str(self.id)+')'

    def __len__(self):
        """Length is determined by the number of dimensions this halfspace is in"""
        return len(self.norm)

# Properties for formatting/sanitizing inputs
    @property
    def rhs(self):
        """RHS of the halfspace"""
        return self._rhs

    @property
    def point(self):
        """Get the defining point of the halfspace"""
        return self._point
    @point.setter
    def point(self,value):
        """Set the defining point of the halfspace. Re-adjusts RHS when setting"""
        self._point=value
        self._rhs=np.dot(self.norm,self._point)

    @property
    def norm(self):
        """Get the normal unit vector"""
        return self._norm
    @norm.setter
    def norm(self,value):
        """Set the normal unit vector by normalizing any given vector"""
        self._norm=np.array(value)/np.linalg.norm(value)

    @property
    def eps(self):
        """Get the detection/error limit"""
        return self._eps
    @eps.setter
    def eps(self,value):
        """Set the detection/error limit"""
        self._eps=min(value,1)
        self._dec=int(max(0,-np.log10(self._eps)))

    @property
    def dec(self):
        """Get the detection/error limit power"""
        return self._dec
    @dec.setter
    def dec(self,value):
        """Set the detection/error limit power"""
        self._dec=max(0,int(value))
        self._eps=10**(-self._dec)

# TODO: Optimize Static Functions (e.g. key) to not require calculation every time