
from .LatticeGraph import LatticeGraph
from .Node import Node
from .Halfspace import Halfspace
from .Search import Search

import traceback
import logging
__log=logging.getLogger('fea.core')

"""
This module contains the flux_envelople_analysis caller function, which will generate a lattice graph of reduced dimensions for a linear optimization problem.
"""

def flux_envelope_analysis(model,variables,max_value=1000,max_iter=1000,eps=10**-4):
    """Run Flux Envelope Analysis [1]_ on a model for the given variables

    Parameters
    ----------
        model : Optlang.Model,
            The original linear program to be reduced
        variables : iterable
            A list of target variables contained in the model
        max_value : positive number
            Maximum/Minimum Value for each variable (-max_value<=variable<=max_value). Will be applied to all variables with bounds greater than limit. Default 1000.
        max_iter : positive integer
            Maximum number of optimization step iterations
        eps : float
            Detection limit. Default 1E-4

    Returns
    -------
        solution : fea.LatticeGraph
            The FEA solution

    Notes
    ------
    This routine will attempt to find a complete solution, but does not guarantee that
    the solution returned will be complete. Always check the 'complete' attribute of the
    solution before utilizing.

    """
    obj=LatticeGraph(model,variables,max_value=max_value,eps=eps)
    res=obj.solve(max_iter)
    return obj

__exports__=[flux_envelope_analysis]