#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 10:19:13 2024

@author: rubenmitchell
"""

#%% import libraries
from docplex.mp.model import Model
#%% import scripts



#%% Model class

class model: 
    def __init__(self, name, data):
        self.data = data
        self.model_instance = Model(name)
        self.x = self.model_instance.binary_var_dict(self.data.A, name = 'x')
        self.model_instance.minimize(self.model_instance.sum(self.x[i,j]*self.data.dists[i,j] for (i,j) in self.data.A))
        self.model_instance.add_constraints(self.model_instance.sum(self.x[p] for p in self.data.A if i in p) == 2 for i in self.data.V)
        self.solution = []
    
    def solve(self):
        self.solution = self.model_instance.solve(log_output = True)