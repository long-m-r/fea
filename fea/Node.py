#!/usr/bin/env python3
import logging
log=logging.getLogger('fea.node')

import itertools

import numpy as np
from .util import lstsq

class Node(frozenset):
    _id_gen = iter(itertools.count())

    def __new__(cls,*args,**kwargs):
        return super().__new__(cls,*args)

    def __init__(self,*args,**kwargs):
        """A container for a single Node in the LatticeGraph as defined in :class:`fea.LatticeGraph`

        Each Node is an immutable set (inherits from :class:`frozenset`) containing solely objects of :class:`fea.Halfspace`.
        The nodes are thus defined by the set of facets which define them.

        Parameters
        ----------
        iterable of :class:`fea.Halfspace` objects
        n : integer
            If defined, specifies the number of dimensions of the problem (necessary only for an empty set where this cannot be found by the facets)
        eps : float
            If defined, the margin of error for this node
        """
        self._id=-1
        self._n=kwargs.pop('n',None)
        self._eps=kwargs.pop('eps',None)
        self._graph=None
        super().__init__()

    def _graph_add(self,graph):
        """internal method to add a graph"""
        self._graph=graph

    def _graph_remove(self):
        """internal method to unset the graph"""
        self._graph=None

    @property
    def complete(self):
        """Boolean indicating whether the node is completed"""
        if self._graph is not None:
            return self._graph.node[self].get('complete',False)
        return False


    @property
    def id(self):
        """Nodes unique ID number. Lazily assigned"""
        # Note that ID's are unique, but assigned lazily. They do not correspond to creation order but to id request order!
        if self._id < 0:
            self._id=next(Node._id_gen)
        return self._id

    @property
    def n(self):
        """The max norm of of the halfspace, indicates level"""
        if self._n is None:
            self._n = max([len(h) for h in self])
        return self._n

    @property
    def eps(self):
        """The detection limit of the method"""
        if self._eps is None:
            self._eps = max([h.eps for h in self])
        return self._eps

    @property
    def real_count(self):
        """The number of real halfspaces in this node"""
        try:
            return self._real
        except AttributeError:
            self._real=np.count_nonzero([h.real for h in self])
            return self._real

    @property
    def real(self):
        """Whether this node is real or not"""
        return self.real_count==len(self)

    @property
    def level(self):
        """The level of this node as an int"""
        try:
            return self._level
        except AttributeError:
            self._level = max([0,self.n-len(self)])
            return self._level

    @property
    def point(self):
        """
        An array of points corresponding to this node as a :class:`numpy.ndarray`
        """
        try:
            return self._point
        except AttributeError:
            if self.level>0 or self.score<0:
                # Can't have a true vertice if we have psuedo facets or aren't a vertex
                if self._graph is not None:
                    return [v.point for v in self._graph.get_vertices(node=self,real=self.real,complete=self.complete)]
                else:
                    return []
            else:
                log.info('Finding Common Point for '+repr(self))
                lst=[h for h in self if h.real]
                A=[h.norm for h in lst]
                b=[h.rhs for h in lst]
                self._point=lstsq(A,b,self.eps)
            return self._point


    @property
    def required_halfspaces(self):
        """Number of halfspaces needed to be complete"""
        req=set()
        for f in self:
            req|=f.required_halfspaces
        return req

    @property
    def valid_domain(self):
        """Whether the domain is valid or not"""
        return len(self.required_halfspaces - self)==0

    def facet_has_vertex(self,other):
        """
        Indicates whether this node contain a vertex

        Parameters
        ----------
        other : :class:`fea.Node`
            The node representing a vertex to check

        Returns
        -------
        True if the vertex is contained in this Node
        """
        if other.point is None:
            raise TypeError(repr(other)+' is not a vertex!')

        for f in self-other:
            if not f.contains(other.point):
                return False

        return True

    def vertex_has_facet(self,other):
        """
        Indicates whether this vertex is contained in another facet

        Parameters
        ----------
        other : :class:`fea.Node`
            The node representing a facet to check

        Returns
        -------
        True if the facet is contained by this vertex
        """
        if self.point is None:
            raise TypeError(repr(self)+' is not a vertex!')

        for f in other-self:
            if not f.contains(self.point):
                return False

        return True

    @property
    def score(self):
        """Adjusted level to show whether the node is complete. Used for sorting"""
        try:
            return self._score
        except AttributeError:
            self._score = self.real_count-(self.n-self.level)
            return self._score

    @property
    def sort_key(self):
        """A tuple indicating how promising this node is for further searches

        We want to search the highest level and lowest scoring nodes in general.
        Lower levels are more specific and provide less overall information.
        Higher scores indicate that the node is further away from being complete
        so the next solution likely will not return an obvious facet.
        """
        # Want best target to have minimum value.
        return (self.level,-self.score)

# Objective Functions
    def orthogonal_vector(self,children=None):
        """Calculate an orthoganol vector to the current node
        
        Parameters
        -----------
        children : iterable of :class:`fea.Halfspace` objects
            If specified, children is an iterable containing facets that we should NOT search towards
        """

        # Find all the child bounding facets from this node (we don't want to find them again after this search!)
        if children is None or len(children)==0:
            children=[self.random_vector()]
        else:
            children=[f.norm for f in children]

        # Figure out the maximum number of children we can incorporate in the solve
        nchild=min(len(children),self.n-len(self))

        # TODO: One shot solve for orthoganol vector! This would be more robust than the current approach. Requires optimization
        # We need to be orthoganol to each facet in this node and non-zero in the direction specified
        a=np.array([h.norm for h in self]+children[0:nchild])
        b=np.array([0]*len(self)+[1]*nchild)

        # Solve for the vector
        if log.getEffectiveLevel()<=10:
            log.info('Solving for orthoganol vector for '+repr(self))
        r=lstsq(a,b,self.eps)

        # Make sure it is not pointing towards any of the children we didn't previously utilize
        for c in children[nchild:]:
            if np.dot(c,r)<0:
                log.error('Could not find an orthogonal vector that points away from all children for '+repr(self))
                raise ValueError('Could not find an orthogonal vector that did not point towards a child.')

        return r/np.linalg.norm(r)

    def random_vector(self):
        """Generate a random search vector
        """
        r=np.random.rand(self.n)-.5
        return r/np.linalg.norm(r)

    # Cast super to node for frozenset operators
    def __and__(self,other):
        return Node(super().__and__(other),n=self.n,eps=self._eps)
    def __or__(self,other):
        return Node(super().__or__(other),n=self.n,eps=self._eps)
    def __xor__(self,other):
        return Node(super().__xor__(other),n=self.n,eps=self._eps)
    def __sub__(self,other):
        return Node(super().__sub__(other),n=self.n,eps=self._eps)

    # Utility Functions
    def __str__(self):
        try:
            return self._str
        except AttributeError:
            try:
                if self.level!=0:
                    raise AttributeError
                self._str=repr(self)+' {level='+str(self.level)+'; score='+str(self.score)+'; point='+str(self._point)+'}'
            except AttributeError:
                self._str=repr(self)+' {level='+str(self.level)+'; score='+str(self.score)+'}'
            for h in self:
                self._str+='\n\t'+str(h)
            return self._str


    def __repr__(self):
        return 'Node('+str(self.id)+')'