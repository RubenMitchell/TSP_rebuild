#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 10:19:13 2024

@author: rubenmitchell
"""

#%% import libraries
from docplex.mp.model import Model
#%% import scripts
from heuristic_functions import heuristics


#%% Model class

class model: 
    def __init__(self, name, data):
        self.data = data
        self.model_instance = Model(name)
        self.x = self.model_instance.binary_var_dict(self.data.A, name = 'x')
        self.model_instance.minimize(self.model_instance.sum(self.x[i,j]*self.data.dists[i,j] for (i,j) in self.data.A))
        self.model_instance.add_constraints(self.model_instance.sum(self.x[p] for p in self.data.A if i in p) == 2 for i in self.data.V)
        self.solution = []
        self.model_instance.parameters.preprocessing.presolve = 0
    
    def solve(self):
        self.solution = self.model_instance.solve(log_output = True)
        
    def build_warmstart(self):
        heu = heuristics(self)
        heu.greedy_blank_sart()
        heu.two_opt()
        warmstart = self.model_instance.new_solution()
        for e in self.data.A:
            value = 0
            if e in heu.current_tour:
                value = 1
            warmstart.add_var_value(self.x[e], value)
        self.model_instance.add_mip_start(warmstart)