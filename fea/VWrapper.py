#!/usr/bin/env python3

class VWrapper:
    _coeff=[1,-1]

    def __init__(self,var,model):
        self.complex=False
        if isinstance(var,VWrapper):
            self.complex=var.complex
            self.vars=[model.variables[v.name] for v in var.vars]
        else:
            try:
                self.vars=[model.variables[var.forward_variable.name],model.variables[var.reverse_variable.name]]
                self.complex=True
            except AttributeError:
                self.vars=[model.variables[var.name]]

    @property
    def expr(self):
        if self.complex:
            return self.vars[0]-self.vars[1]
        return self.vars[0]
        
    @property
    def primal(self):
        return sum([self._coeff[i]*self.vars[i].primal for i in range(len(self.vars))])

    @property
    def dual(self):
        return sum([self._coeff[i]*self.vars[i].dual for i in range(len(self.vars))])

    @property
    def name(self):
        return self.vars[0].name

    @property
    def ub(self):
        if self.complex and self.vars[0].ub==0:
            return -self.vars[1].lb
        return self.vars[0].ub
    @ub.setter
    def ub(self,val):
        if self.complex:
            if val<0:
                self.vars[0].ub=0
                self.vars[1].lb=-val
            else:
                self.vars[0].ub=val
                self.vars[1].lb=0
        else:
            self.vars[0].ub=val

    @property
    def lb(self):
        if self.complex:
            if self.vars[1].ub==0:
                return self.vars[0].lb
            else:
                return -self.vars[1].ub
        else:
            return self.vars[0].lb
    
    @lb.setter
    def lb(self,val):
        if self.complex:
            if val>0:
                self.vars[1].ub=0
                self.vars[0].lb=val
            else:
                self.vars[1].ub=-val
                self.vars[0].lb=0
        else:
            self.vars[0].lb=val

    def __str__(self):
        return str(self.name)