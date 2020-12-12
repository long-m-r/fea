#!/usr/bin/env python3
import logging
log=logging.getLogger('fea.graph')

import numpy as np
import pandas as pd

from networkx import DiGraph, NetworkXError

from itertools import combinations, count, chain, accumulate
from functools import reduce
from sortedcontainers import SortedListWithKey

from .util import lstsq
from .Node import Node
from .Search import Search
from .VWrapper import VWrapper


class LatticeGraph(DiGraph):
    """A directed graph containing the Face Lattice Graph for the reduced space solution
    
    Inherits from :class:`networkx.DiGraph`

    Parameters
    -----------
    problem : :class:`OptLang.Model`
    variables : iterable
        A list of target variables contained in the model
    max_value : positive number
        Maximum/Minimum Value for each variable (-max_value<=variable<=max_value). Will be applied to all variables with bounds greater than limit. Default 1000.
    eps : float
        Detection limit. Default 1E-4
    """
    def __init__(self,problem,variables,max_value=1000,eps=10**-6):
        self.EPS=eps
        self._n=len(variables)

        # Minimum f-vector from n-simplex
        self._minimum_f_vector=np.fromiter(accumulate(range(1,self.N+3), func=lambda p,k: p*(self.N+2-k+1)/(k-1)),dtype=np.int)

        # Get the problem (should work for OptLang, Cameo, and CobraPy>=0.6.0)
        try:
            # This approach should work for Cameo and the latest CobraPy
            self._problem=problem.solver.interface.Model.clone(problem.solver)
        except AttributeError:
            try:
                # This should work for older OptLang model
                self._problem=problem.interface.Model.clone(problem)
            except AttributeError:
                # This should work for newer
                self._problem=problem.clone(problem)

        # Get the variables and set the maximum bounds for search.
        self._variables=[VWrapper(v,self._problem) for v in variables]

        for v in self._variables:
            if v.lb is None or v.lb < -max_value:
                v.lb=-max_value
            if v.ub is None or v.ub > max_value:
                v.ub=max_value

        # TODO: Multi-Threaded Pool of Searchers?
        self.searcher=Search(self._problem,self._variables,eps=self.EPS)
        super().__init__()
        self.reset()

    # def hot_restart(self):
    #         if len(self.queue)>0:
    #             raise ExecutionError('Already have searchers')
            
    #         nodes=self.nodes
    #         for n in nodes:
    #             if n in self:
    #                 if not n.real:
    #                     self.remove_node(n)
    #                 elif not self.node[n].get('complete',False):
    #                     self.queue.add(n)

    def reset(self):
        """Completely clear the graph to its initial state"""
        self.clear()
        self._iterations=0
        self._trace_iter=iter(count())
        self._f_vector=[0 for i in range(self.N+1)]
        self._complete_halfspaces=set()
        self.queue=SortedListWithKey([],key=lambda n: n.sort_key)

        # Seed the graph with the Polytope Node (empty set)
        self._polytope_node=Node(n=self.N,eps=self.EPS)
        self.add_node(self._polytope_node, trace=0)

        
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Basic Graph Properties <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    @property
    def M(self):
        """Original problem dimensions"""
        return len(self._problem.variables)  # Note, this may be twice as large if variables are in standard form (i.e. always positive)
    @property
    def N(self):
        """Reduced problem/graph dimensions"""
        return self._n
    @property
    def POLYTOPE_LEVEL(self):
        """Level of the Polytope Node"""
        return self.N
    @property
    def FACET_LEVEL(self):
        """Level of the Facet Nodes"""
        return self.N-1
    @property
    def EDGE_LEVEL(self):
        """Level of the Edge Nodes (*i.e.* 1)"""
        return 1
    @property
    def VERTEX_LEVEL(self):
        """Level of the Vertex Nodes (*i.e.* 0)"""
        return 0
    @property
    def _trace(self):
        return next(self._trace_iter)
    
    @property
    def EPS(self):
        """The detection limit of the method"""
        return self._eps
    @property
    def DEC(self):
        """The order-of-magnitude of EPS"""
        return self._dec
    
    @EPS.setter
    def EPS(self,value):
        self._eps=min(value,1)
        self._dec=int(max(0,-np.log10(self._eps)))
    @DEC.setter
    def DEC(self,value):
        self._dec=max(0,int(value))
        self._eps=10**(-self._dec)

    @property
    def f_vector(self):
        """The f-vector for the graph

        Notes
        -----
        This vector is continously updated and maintained and is NOT re-calculated
        """
        return np.array([1]+self._f_vector, dtype=np.int)

    @property
    def modified_euler_characteristic(self):
        """The modified Euler characteristic

        Notes
        -----
        A completed and enclosed polytope will always have a modified Euler characteristic of zero. The result is equivalent to:

        .. math::
            (-1)^{N} * \sum^{N}_{i=-1} (-1)^{i} f_{i}

        """
        return reduce(lambda x,y: -x+y, self.f_vector)

    @property
    def complete(self):
        return np.all(np.greater_equal(self.f_vector,self._minimum_f_vector)) and self.modified_euler_characteristic==0
        
    # def __str__(self):
    #     nodes=sorted(self.nodes(),key=lambda x: x.sort_key, reverse=True)
    #     rstr='\n'.join([str(n).replace('\n','\n\t') for n in nodes])
    #     return 'Lattice Graph:\n'+rstr

#>>>>>>>>>>>>>>>>>>>>>>>> Graph Modification Functions <<<<<<<<<<<<<<<<<<<<<<<<
    def add_node(self, node, **kwargs):
        """Add a node to the graph

        Parameters
        ----------
        node : :class:`fea.Node`
            The node to add. It must be a Node object.

        Notes
        -----
        Other named arguments may be passed as attributes of the node in the graph (see :func:`networkx.DiGraph.add_node()`)

        All added nodes will be automatically connected to appropriate parent and child nodes. The f-vector will be updated accordingly
        """
        if not isinstance(node,Node) or node in self or not node.valid_domain:
            raise ValueError('Invalid node for addition:\n'+str(node))
        if 'trace' not in kwargs:
            kwargs['trace']=self._trace
        is_recurse=kwargs.pop('_recurse',False)

        # Additional checks required for adding vertices
        if node.level==self.VERTEX_LEVEL and node.real:
            other_facets=set()
            for f in self.get_facets(complete=None):
                if node.vertex_has_facet(f):
                    other_facets|=f
                other_facets-=node
            if len(other_facets)>0:
                node = node | other_facets
                log.info('Adding facets '+','.join([repr(f) for f in other_facets])+' to form '+repr(node))
                
                for v in self.get_vertices(complete=None):
                    if v<=node:
                        log.info('Removing subset vertex '+repr(v)+' in favor of '+repr(node))
                        self.remove_node(v)


        log.info('Adding '+repr(node))
        node._graph_add(self)
        super().add_node(node, **kwargs)
        
        if node.level>self.VERTEX_LEVEL:
            self.queue.add(node)

        # Search for all possible parents and connect to them. Be very careful with psuedo facets. We *MUST NOT* propogate psuedo facets UP the graph!
        if node.level<self.N:
            for potential_parent in combinations(node, self.N-node.level-1):
                pnode=Node(potential_parent, n=self.N, eps=self.EPS)
                if pnode in self:
                    if not self.has_edge(pnode,node):
                        self.add_edge(pnode, node, trace=kwargs['trace'])
                elif pnode.valid_domain:
                    try:
                        self.add_node(pnode, trace=kwargs['trace'], _recurse=True)
                        self.add_edge(pnode, node, trace=kwargs['trace'])
                    except ValueError:
                        # The node was invalid for some reason, don't really care why
                        pass

        # If it's a real facet, then check all vertices and add accordingly
        if node.level==self.FACET_LEVEL and node.real:
            facet, = node
            for v in self.get_vertices(complete=None):
                if not node <= v and node.facet_has_vertex(v):
                    log.info('Updating vertex '+repr(v)+' to include '+repr(facet))
                    nodedict = self.node[v]
                    nodedict['_recurse']=True
                    self.remove_node(v, _recurse=True)
                    self.add_node(node | v, **nodedict)

        # Now that we're done making changes, check that graph completeness hasn't changed
        if not is_recurse:
            self._update_graph_completeness()
            self._update_node_completeness(node)

        return node

    def remove_node(self, node, **kwargs):
        """Remove a node from the graph

        Parameters
        ----------
        node : :class:`fea.Node`
            The node to remove. It must look like a Node object.

        Notes
        -----
        All added nodes will be automatically disconnected and orphaned child nodes will be removed automatically. The f-vector will be updated accordingly
        """
        is_recurse=kwargs.get('_recurse',False)
        node._graph_remove()

        if not is_recurse:
            node=self.get_node(node)

        # Remove all successor nodes dependent upon this one BEFORE we delete this node
        for n in self.successors(node):
            try:
                self.remove_node(n, _recurse=True)
            except ValueError:
                pass # It's already been deleted

        # Remove this node
        log.info('Removing '+repr(node))
        if node not in self:
            raise ValueError('Invalid node for deletion:\n'+str(node))
        elif self.node[node].get('complete',False):
            self._f_vector[node.level]-=1
            pred_inherit=self.predecessors(node,real=True)

            super().remove_node(node)
            self.queue.discard(node)
            
            for p in pred_inherit:
                self.node[p]['_complete_children'] = min(0,self.node[p].get('_complete_children',0)-1)
                self._update_node_completeness(p)

            if not is_recurse:
                self._update_graph_completeness()
        else:
            super().remove_node(node)
            self.queue.discard(node)
            
            if not is_recurse:
                self._update_graph_completeness()
        
        # Now that we're done making changes, check that graph completeness hasn't changed
        if not is_recurse:
            self._update_graph_completeness()

    def add_edge(self,nodefrom,nodeto,**kwargs):
        """Connect a parent node (higher level) to a child node (lower level).
        
        Parameters
        ----------
        nodefrom, nodeto : :class:`fea.Node`
            Nodes to connect in lattice graph
        """
        if log.getEffectiveLevel() <= 10:
            # Only need the actual nodes if we're logging. Kind of expensive hence the check
            nodefrom, nodeto = self.get_node(nodefrom), self.get_node(nodeto)
            log.info('Adding edge '+repr(nodefrom)+'->'+repr(nodeto))

        if not kwargs.get('trace',False):
            kwargs['trace']=self._trace
        super().add_edge(nodefrom, nodeto, **kwargs)

        if self.node[nodeto].get('complete',False):
            self.node[nodefrom]['_complete_children']=self.node[nodefrom].get('_complete_children',0)+1
            self._update_node_completeness(nodefrom)

#>>>>>>>>>>>>>>>>>> Graph Query Functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    def successors(self,node,real=None,complete=None):
        return list(filter(lambda n: (real is None or n.real==real) and (complete is None or self.node[n].get('complete',False)==complete), super().successors(node)))

    def predecessors(self,node,real=None,complete=None):
        return list(filter(lambda n: (real is None or n.real==real) and (complete is None or self.node[n].get('complete',False)==complete), super().predecessors(node)))

    def get_node(self, node):
        if isinstance(node,int):
            for n in self.nodes_iter():
                if n.id==node:
                    return n
        elif node in self:
            # Requires iterating the whole stinkin graph in the worst case
            for n in self.node.keys():
                if n==node:
                    return n
        return None

    def get_nodes_of_level(self,level,connected=None,real=None,complete=None):
        res=set()

        if connected is not None:
            if connected.level<level:
                for n in self.predecessors(connected,real=real,complete=complete):
                    if n.level==level:
                        res.add(n)
                    elif n.level<level:
                        res|=self.get_nodes_of_level(level,connected=n,real=real,complete=complete)
            elif connected.level>level:
                for n in self.successors(connected,real=real,complete=complete):
                    if n.level==level:
                        res.add(n)
                    elif n.level>level:
                        res|=self.get_nodes_of_level(level,connected=n,real=real,complete=complete)
            elif connected.level==level:
                if (real is None or connected.real==real) and (complete is None or self.node[connected].get('complete',False)==complete):
                    res.add(connected)
        else:
            for n in self.nodes_iter():
                if n.level==level and (real is None or n.real==real) and (complete is None or self.node[n].get('complete',False)==complete):
                    res.add(n)
        return res

    def get_vertices(self,node=None,real=True,complete=True,pandas=False):
        # TODO: Maintain list of these? This is an expensive step that is oft repeated
        if pandas:
            res=[]
            for n in self.get_nodes_of_level(self.VERTEX_LEVEL,connected=node,real=real,complete=complete):
                res+=[[n.id]+list(np.round(n.point,self.DEC))]
            return pd.DataFrame(data=res,columns=['id']+[v.name for v in self._variables])
        return self.get_nodes_of_level(self.VERTEX_LEVEL,connected=node,real=real,complete=complete)

    def get_facets(self,node=None,real=True,complete=True,pandas=False,rhs_label='RHS'):
        # TODO: Maintain list of these? This is an expensive step that is oft repeated
        if pandas:
            res=[]
            for n in self.get_nodes_of_level(self.FACET_LEVEL,connected=node,real=real,complete=complete):
                hp=next(iter(n))
                res+=[[hp.id]+list(np.round(np.append(hp.norm,hp.rhs),self.DEC))]
            return pd.DataFrame(data=res,columns=['_id']+[v.name for v in self._variables]+[rhs_label])
        return self.get_nodes_of_level(self.FACET_LEVEL,connected=node,real=real,complete=complete)

#>>>>>>>>>>>>>>>>> Functions for Searching/Solving the Graph <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    def solve(self,max_iter=50,exhaust=False):
        ind=0
        while len(self.queue)>0 and ind<max_iter and (exhaust or (self.queue[0].level==1 and self.queue[0].real) or not self.complete):
            newsearch=self.queue[0]
            res=self.search(newsearch)

            if res:
                log.info('Done @ '+str(ind)+'\n\n')
                ind+=1
            else:
                log.info('Done & Exhausted Node @ '+str(ind)+'\n\n')
                self.queue.discard(newsearch)
                
        if log.getEffectiveLevel()<=10:
            log.info('Complete:'+str(self.complete)+' QueueLength:'+str(len(self.queue)))
            log.info('Finished Search after '+str(ind)+' iterations with f_vector:'+str(self.f_vector)+'='+str(self.modified_euler_characteristic))
            log.info(str(self))

        self._iterations+=ind
        return ind

    def search(self,node=None):
        if node is None:
            node=self.queue[0]

        if log.getEffectiveLevel()<=10:
            log.info('Searching { f_vector:'+str(self.f_vector)+'='+str(self.modified_euler_characteristic)+'} for '+repr(node))

        # If it's a vertex
        if node.level==0:
            log.info('Vertex. Aborting Solve.')
            return False
        elif node.level==1 and len([s for s in self.successors(node,real=True)])>=2:
            log.info('Completed Edge. Aborting Solve.')
            return False

        # First, find any facets that we already know exist for this node
        knownfacets=set()
        for s in self.successors(node):
            if s.real or self.edge[node][s].get('searched',False):
                knownfacets|=s
            else:
                knownfacets|=set([f for f in s if f.real])
        knownfacets-=node
        log.info('Found bounding '+','.join([repr(k) for k in knownfacets])+' for '+repr(node))


        # Now, look for an objective direction
        try:
            obj=node.orthogonal_vector(knownfacets)
        except ValueError:
            # We couldn't find a valid/new objective direction. We're done with this node
            log.error('Could not find an orthogonal vector to child nodes!')
            return False

        # Now, we can solve the Problem
        trace=self._trace
        log.info('Solving with trace '+str(trace))
        self.searcher.set(obj,node)
        optimal=self.searcher.get_solution()
        if not optimal:
            log.error('Solver Error. Aborting.')
            return False

        # Calculate the new bounding halfspace and make sure we have the real one
        halfspace = self.searcher.bounding_halfspace()
        log.info('Found '+str(halfspace))
        facet=Node([halfspace])
        if facet in self:
            halfspace, = self.get_node(facet)
        log.info('Changing to '+str(halfspace))

        child_node=self.add_node( Node( node | set([halfspace]) ), trace=trace)
        self.edge[node][child_node]['searched']=trace

        return True

#>>>>>>>>>>>>>>>>>>>> Graph and Node Completion Functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
    def _update_graph_completeness(self,**kwargs):
        iteration = kwargs.get('_iteration',0)
        update=False

        if iteration==0:
            complete_halfspaces, incomplete_halfspaces = set(), set(chain.from_iterable(self.get_facets(real=True,complete=None)))
            complete_vertices, incomplete_vertices = set(self.get_vertices(real=True,complete=None)) , set()

            # Check Potential Halfspaces for Enough Vertices
            for ih in incomplete_halfspaces:
                # Get number of vertices that include it
                if len([1 for cv in complete_vertices if ih in cv])>=self.N:
                    complete_halfspaces.add(ih)
            incomplete_halfspaces-=complete_halfspaces
            
        else:
            complete_halfspaces, incomplete_halfspaces = kwargs['_halfspaces']
            complete_vertices, incomplete_vertices = kwargs['_vertices']

        # Check Potential Vertices for Enough Halfspaces
        for cv in complete_vertices:
            if len(cv & complete_halfspaces) < self.N:
                incomplete_vertices.add(cv)
                update=True
        complete_vertices-=incomplete_vertices
        
        # Check Complete Halfspaces for Enough Vertices
        for ch in complete_halfspaces:
            # Get number of vertices that include it
            if len([1 for cv in complete_vertices if ch in cv])<self.N:
                incomplete_halfspaces.add(ch)
                update=True
        complete_halfspaces-=incomplete_halfspaces

        # log.info('Checking completeness iteration '+str(iteration)+'\n'+str())
        
        if update and len(complete_halfspaces)>0 and len(complete_vertices)>0:
            self._update_graph_completeness(_halfspaces=(complete_halfspaces,incomplete_halfspaces),_vertices=(complete_vertices,incomplete_vertices),_iteration=iteration+1)
        else:
            self._complete_halfspaces=complete_halfspaces
            for v in chain(complete_vertices,incomplete_vertices):
                self._update_node_completeness(v)
            
    def _update_node_completeness(self, node):
        current = self.node[node].get('complete',False)
        possible = node.real and len(node & self._complete_halfspaces) >= (self.N - node.level) and (node.level==self.VERTEX_LEVEL or self.node[node].get('_complete_children',0)>node.level)
        
        if current != possible:
            #Children may be screwed up, make sure they're ok            
            if possible:
                self.node[node]['complete']=True
                self._f_vector[node.level]+=1
                for p in self.predecessors(node,real=True):
                    self.node[p]['_complete_children']=self.node[p].get('_complete_children',0)+1
                    self._update_node_completeness(p)
            else:
                self.node[node]['complete']=False
                self._f_vector[node.level]-=1
                for p in self.predecessors(node,real=True):
                    self.node[p]['_complete_children']=min(0,self.node[p].get('_complete_children',0)-1)
                    self._update_node_completeness(p)

#>>>>>>>>>>>>>>>>>>>>>>> Output/Formatting Functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    def to_optlang_model(self,replace_variables=True):
        interface=self._problem.interface
        newmodel=interface.Model()

        if replace_variables:
            variables=[interface.Variable(v.name,ub=v.ub,lb=v.lb) for v in self._variables]
            for h in self._complete_halfspaces:
                newmodel.add(interface.Constraint(np.dot(variables,h.norm),lb=h.rhs,ub=None,name=repr(h)))
        else:
            for h in self._complete_halfspaces:
                new=h._ol_new_constraint(newmodel,eps=0)
                new.ub=None
                new.name=repr(h)
                newmodel.add(new)

        return newmodel



    def to_scipy_optimize_model(self):
        kwargs={}


        def fun(x,h):
            return np.dot(x,h.norm)-h.rhs
        def jac(x,h):
            return h.norm

        # Populate constraints
        constraints=[]
        for fn in self.get_facets(real=True):
            constraints.append({
                'type':'ineq',
                'fun': fun,
                'jac': jac,
                'args': [next(iter(fn))]})
        kwargs['constraints']=constraints

        # Populate other optional/required args
        kwargs['x0']=self.get_vertices(real=True)[0].point
        kwargs['bounds']=[(v.lb,v.ub) for v in self._variables]
        kwargs['tol']=self.EPS

        return kwargs